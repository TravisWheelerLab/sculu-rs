use anyhow::{anyhow, bail, Result};
use assert_cmd::Command;
use chrono::Duration;
use clap::{builder::PossibleValue, Parser, ValueEnum};
use csv::{ReaderBuilder, WriterBuilder};
use itertools::Itertools;
use kseq::parse_reader;
use log::{debug, info};
use newick::Newick;
use rand::{seq::SliceRandom, thread_rng};
use regex::Regex;
use serde::{Deserialize, Serialize};
use std::{
    cmp,
    collections::{HashMap, HashSet},
    ffi::OsStr,
    fmt,
    fs::{self, File},
    io::{BufRead, BufReader, BufWriter, Write},
    path::{Path, PathBuf},
    time::Instant,
};
use walkdir::WalkDir;
use which::which;

//#[derive(Parser, Debug)]
//#[command(arg_required_else_help = true, version, about)]
//pub struct Cli {
//    #[command(subcommand)]
//    pub command: Option<Command>,
//}
//
//#[derive(Parser, Debug)]
//pub enum Command {
//    /// Create sufr file
//    Cluster(ClusterArgs),
//
//    Run(RunArgs)
//}

/// SCULU subfamily clustering tool
#[derive(Debug, Parser)]
#[command(author, version, about)]
pub struct Args {
    /// FASTA file of subfamily consensi
    #[arg(long, value_name = "CONSENSI")]
    pub consensi: PathBuf,

    /// Directory of instance files for each subfamily
    #[arg(long, value_name = "INSTANCES", required = true, num_args=1..)]
    pub instances: Vec<PathBuf>,

    /// Output directory
    #[arg(short, long, value_name = "OUTDIR", default_value = "sculu-out")]
    pub outdir: PathBuf,

    /// Lambda value
    #[arg(long, value_name = "LAMBDA", default_value = "0.1227")]
    pub lambda: f64,

    /// Independence threshold
    #[arg(short, long, value_name = "IND", default_value = "0.5")]
    pub independence_threshold: f64,

    /// Confidence margin
    #[arg(short, long, value_name = "CONF", default_value = "3")]
    pub confidence_margin: isize,

    /// Path to rmblastn
    #[arg(long, value_name = "ALIGNER")]
    pub aligner: Option<String>,

    /// Path to RepeatModeler/Refiner
    #[arg(long, value_name = "REFINER")]
    pub refiner: Option<String>,

    /// Number of threads for rmblastn/Refiner
    #[arg(long, value_name = "THREADS")]
    pub num_threads: Option<usize>,

    /// PERL5LIB, e.g., to find RepeatMasker/RepeatModeler
    #[arg(long, value_name = "PERL5LIB")]
    pub perl5lib: Option<String>,

    /// Path to rmblastn
    #[arg(long, value_name = "RMBLAST_DIR")]
    pub rmblast_dir: Option<String>,

    /// Alignment matrix
    #[arg(long, value_name = "MATRIX")]
    pub alignment_matrix: Option<PathBuf>,

    /// Alignment gap init
    #[arg(long, value_name = "GAPINIT", default_value = "-25")]
    pub align_gap_init: i64,

    /// Alignment extension
    #[arg(long, value_name = "EXT", default_value = "-5")]
    pub align_extension: i64,

    /// Alignment minimum match
    #[arg(long, value_name = "MINMATCH", default_value = "7")]
    pub align_min_match: i64,

    /// Alignment bandwidth
    #[arg(long, value_name = "BANDWIDTH", default_value = "14")]
    pub align_bandwidth: i64,

    /// Alignment mask level
    #[arg(long, value_name = "MASKLEVEL", default_value = "101")]
    pub align_mask_level: i64,

    /// Alignment minimum score
    #[arg(long, value_name = "MINSCORE", default_value = "200")]
    pub align_min_score: i64,

    /// Log level
    #[arg(short, long)]
    pub log: Option<LogLevel>,

    /// Log file
    #[arg(long)]
    pub log_file: Option<String>,
}

#[derive(Debug, Clone)]
pub enum LogLevel {
    Info,
    Debug,
}

impl ValueEnum for LogLevel {
    fn value_variants<'a>() -> &'a [Self] {
        &[LogLevel::Info, LogLevel::Debug]
    }

    fn to_possible_value<'a>(&self) -> Option<PossibleValue> {
        Some(match self {
            LogLevel::Info => PossibleValue::new("info"),
            LogLevel::Debug => PossibleValue::new("debug"),
        })
    }
}

#[derive(Debug, Deserialize)]
struct RmBlastOutput<'a> {
    score: u32,
    target: &'a str,
    query: &'a str,
}

#[derive(Debug, Serialize, Deserialize)]
struct AlignmentScore {
    score: u32,
    target: String,
    query: String,
}

#[derive(Debug, PartialEq)]
struct Independence {
    f1: String,
    f2: String,
    val: f64,
}

#[derive(Debug, Eq, PartialEq, Hash)]
struct StringPair(String, String);

impl fmt::Display for StringPair {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}/{}", self.0, self.1)
    }
}

#[derive(Debug)]
struct Winners {
    clear_winners: HashMap<String, u32>,
    winning_sets: HashMap<StringPair, u32>,
}

#[derive(Debug)]
struct MergeFamilies<'a> {
    family1: String,
    family2: String,
    outdir: PathBuf,
    instances_100_dir: &'a PathBuf,
    num_threads: usize,
    refiner: &'a Option<String>,
    log: &'a Option<LogLevel>,
    rmblast_dir: &'a Option<String>,
    perl5lib: &'a Option<String>,
}

// --------------------------------------------------
pub fn run(mut args: Args) -> Result<()> {
    let start = Instant::now();

    env_logger::Builder::new()
        .filter_level(match args.log {
            Some(LogLevel::Debug) => log::LevelFilter::Debug,
            Some(LogLevel::Info) => log::LevelFilter::Info,
            _ => log::LevelFilter::Off,
        })
        .target(match args.log_file {
            // Optional log file, default to STDOUT
            Some(ref filename) => env_logger::Target::Pipe(Box::new(
                BufWriter::new(File::create(filename)?),
            )),
            _ => env_logger::Target::Stdout,
        })
        .init();

    debug!("args = {args:#?}");
    if !args.outdir.is_dir() {
        fs::create_dir_all(&args.outdir)?;
    }

    //align_consensi_to_self(&args)?;
    //panic!("stop");

    // Take the longest 100 of the instances
    let instances_100_dir = args.outdir.join("instances_100");
    debug!(
        "Writing 100 longest instances to '{}'",
        instances_100_dir.display()
    );
    fs::create_dir_all(&instances_100_dir)?;
    let instances_100 = take_longest_100(&args, &instances_100_dir)?;
    debug!("instances_100 = '{instances_100:#?}'");

    // Check that the "consensi" sequences have clear relationships
    // to "instances" files. Returns a mutable hashmap from the
    // family names to their instance files for merges.
    debug!("Checking family names/instances");
    let mut family_to_path =
        check_family_instances(&args.consensi, &instances_100)?;
    debug!("family_to_path = '{family_to_path:#?}'");

    // Concatenate the longest 100 into a single file for alignment
    let all_seqs_path = &args.outdir.join("all_seqs.fa");
    debug!(
        "Concatenating all instance sequences into '{}'",
        all_seqs_path.display()
    );

    // How do we select the instances? Use all or a subset?
    //cat_sequences(&instances_100, all_seqs_path)?;

    // Alternative is to use *all* the sequences
    cat_sequences(&args.instances, all_seqs_path)?;

    // Create a starting directory for the merge iterations.
    let top_outdir = args.outdir.clone();
    let round0_dir = top_outdir.join("round000");
    fs::create_dir_all(&round0_dir)?;

    // As the consensi are merged, the family names are concatenated
    // into Newick-formatted strings to track the merges.
    // These can get quite long, causing `makeblastdb` to fail.
    // Make the sequence IDs simple integers by numbering
    // and move the family names to the description.
    let consensi_path = round0_dir.join("consensi.fa");
    debug!("Writing consensi to '{}'", consensi_path.display());
    let mut consensi_seqs = number_consensi(&args.consensi, &consensi_path)?;
    args.consensi = consensi_path;

    // Start merging
    let trailing_semi = Regex::new(";$").unwrap();
    let mut prev_scores: Option<PathBuf> = None;
    for round in 1.. {
        info!(">>> Round {round} <<<");
        args.outdir = top_outdir.join(format!("round{round:03}"));
        fs::create_dir_all(&args.outdir)?;

        // Align all the sequences to the current consensi.
        // On the first round, all the original consensi will be included.
        // On future rounds, only the newly merged consensi will be present.
        let alignment_file = align(&args, all_seqs_path)?;

        // Extract the scores from the alignment file.
        // On the first round, there will be no "prev_scores" file.
        // On the following rounds, the previous round's scores will be
        // included for those families that were not merged.
        let scores_file = extract_scores(
            &alignment_file,
            &prev_scores,
            &args.consensi,
            &args.outdir,
        )?;

        // Save this round's alignment scores for the next
        prev_scores = Some(scores_file.clone());

        // Find the clear/ambiguous winners from the scores
        let winners =
            call_winners(&scores_file, args.lambda, args.confidence_margin)?;

        // Determine independence of all pairs
        let pair_independence = independence(winners);

        // Create a lookup for the scores
        let mut score_lookup: HashMap<StringPair, f64> = HashMap::new();
        for pair in &pair_independence {
            score_lookup.insert(
                StringPair(pair.f1.to_string(), pair.f2.to_string()),
                pair.val,
            );
        }

        // Select only those pair lacking independence
        let non_independent: Vec<_> = pair_independence
            .iter()
            .filter(|v| v.val < args.independence_threshold)
            .collect();

        debug!("{} family pair not independent.", non_independent.len());
        debug!("{non_independent:#?}");

        // Stop the loop when all families are independent
        if non_independent.is_empty() {
            break;
        }

        // Merge the least independent families.
        // The new consensi file will only contain the merged families.
        let new_consensi_path = args.outdir.join("new-consensi.fa");
        let mut new_consensi = open_for_write(&new_consensi_path)?;
        let mut merge_num = 0;
        let mut already_merged: HashSet<String> = HashSet::new();
        for pair in non_independent {
            // The family names may be in Newick format
            let fams1 = parse_newick(&pair.f1);
            let fams2 = parse_newick(&pair.f2);

            // All the family names in this merge
            let mut all_families: Vec<_> = fams1.clone();
            for f2 in &fams2 {
                all_families.push(f2.clone());
            }

            // We might see Fam1->Fam2 and later Fam2->Fam1.
            // Or we might later see Fam1->Fam3.
            // Only merge a family once.
            if all_families.iter().any(|f| already_merged.contains(f)) {
                debug!("Already merged one of {}", all_families.join(", "));
                continue;
            }

            // Increment the number of merged pairs
            merge_num += 1;

            // See if the pair has a score from the other direction.
            // The merges should happen in order from least independent
            // to greater, so this A/B merge "val" will be lower than the
            // symmetrical B/A score.
            let other_score = score_lookup
                .get(&StringPair(pair.f2.to_string(), pair.f1.to_string()))
                .map_or("".to_string(), |v| format!("/{v:0.04}"));

            info!(
                "{}: Merge {} => {} Ind: {:0.04}{other_score}",
                merge_num, &pair.f1, &pair.f2, pair.val
            );

            // Place all merge artefacts into a directory.
            let outdir = args.outdir.join(format!("msa-{merge_num:02}"));
            let new_consensus_seq = merge_families(
                MergeFamilies {
                    family1: pair.f1.clone(),
                    family2: pair.f2.clone(),
                    outdir,
                    instances_100_dir: &instances_100_dir,
                    num_threads: args.num_threads.unwrap_or(num_cpus::get()),
                    refiner: &args.refiner,
                    log: &args.log,
                    rmblast_dir: &args.rmblast_dir,
                    perl5lib: &args.perl5lib,
                },
                &mut family_to_path,
            )?;

            let new_family_newick = format!(
                "({},{}):{:0.02};",
                trailing_semi.replace(&pair.f1, ""),
                trailing_semi.replace(&pair.f2, ""),
                pair.val
            );

            // Write the merged families to the new consensi file
            writeln!(
                new_consensi,
                ">{merge_num} {new_family_newick}\n{new_consensus_seq}",
            )?;

            let f1_len = consensi_seqs.get(&pair.f1).map_or(0, |v| v.len());
            let f2_len = consensi_seqs.get(&pair.f2).map_or(0, |v| v.len());
            let new_len = new_consensus_seq.len();

            info!(
                "{} len was {f1_len}, {} len was {f2_len}, \
                    new consensi len is {new_len}",
                pair.f1, pair.f2
            );

            // Update the consensi mapping
            consensi_seqs.remove(&pair.f1);
            consensi_seqs.remove(&pair.f2);
            consensi_seqs
                .insert(new_family_newick, new_consensus_seq.clone());

            // Keep track of all the families that have been merged
            for family in all_families {
                already_merged.insert(family);
            }
        }

        info!(
            "Wrote {merge_num} consensi to '{}'",
            new_consensi_path.display()
        );
        args.consensi = new_consensi_path;
    }

    let final_consensi_path = args.outdir.join("final-consensi.fa");
    let mut final_consensi = open_for_write(&final_consensi_path)?;
    for (i, (name, seq)) in consensi_seqs.iter().enumerate() {
        writeln!(final_consensi, ">{} {name}\n{seq}", i + 1,)?;
    }

    println!(
        "Finished in {}. Wrote {} final consensi to '{}'",
        format_seconds(start.elapsed().as_secs()),
        consensi_seqs.len(),
        final_consensi_path.display(),
    );

    Ok(())
}

// --------------------------------------------------
fn align_consensi_to_self(args: &Args) -> Result<()> {
    let outfile = args.outdir.join("consensi_cluster");

    let rmblastn_args = vec![
        "-num_alignments".to_string(),
        "9999999".to_string(),
        "-mask_level".to_string(),
        args.align_mask_level.to_string(),
        "-min_raw_gapped_score".to_string(),
        args.align_min_score.to_string(),
        "-num_threads".to_string(),
        args.num_threads.unwrap_or(num_cpus::get()).to_string(),
        "-db".to_string(),
        args.consensi.to_string_lossy().to_string(),
        "-query".to_string(),
        args.consensi.to_string_lossy().to_string(),
        "-out".to_string(),
        outfile.to_string_lossy().to_string(),
        "-outfmt".to_string(),
        "6 score qseqid sseqid".to_string(),
        // TODO: Verify these options
        "-gapopen".to_string(),
        "20".to_string(),
        "-gapextend".to_string(),
        "5".to_string(),
        "-word_size".to_string(),
        "7".to_string(),
        "-xdrop_ungap".to_string(),
        "400".to_string(),
        "-xdrop_gap_final".to_string(),
        "200".to_string(),
        "-xdrop_gap".to_string(),
        "100".to_string(),
        "-dust".to_string(),
        "no".to_string(),
        // TODO: verify remove
        //"-complexity_adjust".to_string(),
    ];

    let rmblastn = match &args.aligner {
        Some(path) => PathBuf::from(path.to_string()),
        _ => which("rmblastn")?,
    };
    let mut cmd = Command::new(&rmblastn);
    if let Some(perl5lib) = &args.perl5lib {
        cmd.env("PERL5LIB", perl5lib);
    }

    // Needed?
    //if let Some(matrix) = &args.alignment_matrix {
    //    let filename = matrix.file_name().unwrap_or_else(|| {
    //        panic!("Failed to get filename from '{}'", matrix.display())
    //    });
    //
    //    let dirname = matrix.parent().unwrap_or_else(|| {
    //        panic!("Failed to get dirname from '{}'", matrix.display())
    //    });
    //
    //    cmd.env("BLASTMAT", dirname);
    //    rmblastn_args.extend_from_slice(&[
    //        "-matrix".to_string(),
    //        filename.to_string_lossy().to_string(),
    //    ]);
    //}

    debug!(
        "Running '{} {}'",
        rmblastn.display(),
        rmblastn_args.join(" ")
    );

    let start = Instant::now();
    let res = cmd.args(rmblastn_args).output()?;
    if !res.status.success() {
        bail!(String::from_utf8(res.stderr)?);
    }

    info!(
        "Consensi self-alignment finished in {}",
        format_seconds(start.elapsed().as_secs())
    );
    Ok(())
}

// --------------------------------------------------
// TODO: remove
//fn _align(args: &Args, all_seqs_path: &Path) -> Result<PathBuf> {
//    // Perform alignment to query
//    let mut align_args = vec![
//        // Don't do complexity adjustment
//        "-raw".to_string(),
//        "-rmblast".to_string(),
//        "-threads".to_string(),
//        args.threads.to_string(),
//        "-gap_init".to_string(),
//        args.align_gap_init.to_string(),
//        "-extension".to_string(),
//        args.align_extension.to_string(),
//        "-minmatch".to_string(),
//        args.align_min_match.to_string(),
//        "-bandwidth".to_string(),
//        args.align_bandwidth.to_string(),
//        "-masklevel".to_string(),
//        args.align_mask_level.to_string(),
//        "-minscore".to_string(),
//        args.align_min_score.to_string(),
//    ];

//    if let Some(matrix) = &args.alignment_matrix {
//        align_args.extend_from_slice(&[
//            "-matrix".to_string(),
//            matrix.to_string_lossy().to_string(),
//        ]);
//    }

//    align_args.extend_from_slice(&[
//        "-alignments".to_string(),
//        all_seqs_path.to_string_lossy().to_string(),
//        args.consensi.to_string_lossy().to_string(),
//    ]);

//    let start = Instant::now();
//    let aligner = match &args.aligner {
//        Some(path) => PathBuf::from(path.to_string()),
//        _ => which("align.pl")?,
//    };

//    debug!("Running '{} {}'", aligner.display(), align_args.join(" "));

//    let mut cmd = Command::new(&aligner);
//    if let Some(perl5lib) = &args.perl5lib {
//        cmd.env("PERL5LIB", perl5lib);
//    }

//    let res = cmd.args(align_args).output()?;
//    if !res.status.success() {
//        bail!(String::from_utf8(res.stderr)?);
//    }

//    info!(
//        "Alignment finished in {}",
//        format_seconds(start.elapsed().as_secs())
//    );

//    let alignment_outfile = args.outdir.join("alignment.txt");
//    let mut output = open_for_write(&alignment_outfile)?;
//    output.write_all(&res.stdout)?;

//    Ok(alignment_outfile)
//}

// --------------------------------------------------
fn bitscore_to_confidence(vals: &[&u32], lambda: f64) -> Result<Vec<f64>> {
    let converted: Vec<_> =
        vals.iter().map(|&&v| (v as f64 * lambda).exp2()).collect();
    let total: f64 = converted.iter().sum();
    if total > 0. {
        Ok(converted.into_iter().map(|v| v / total).collect())
    } else {
        bail!("Sum of converted values equals zero")
    }
}

// --------------------------------------------------
fn align(args: &Args, all_seqs_path: &Path) -> Result<PathBuf> {
    let blast_suffixes: HashSet<_> = [
        "ndb", "nhr", "nin", "njs", "nog", "nos", "not", "nsq", "ntf", "nto",
    ]
    .iter()
    .map(OsStr::new)
    .collect();

    let blast_db_dir = &args.consensi.parent().expect("blast db dir");
    let blast_files: Vec<_> = WalkDir::new(blast_db_dir)
        .max_depth(0)
        .into_iter()
        .filter_map(Result::ok)
        .filter(|e| {
            // Find files
            e.file_type().is_file()
                // With a nonzero size
                && e.metadata().unwrap().len() > 0
                // With a BLAST db suffix
                && blast_suffixes.contains(e.path().extension().unwrap())
        })
        .collect();

    // See if we need to make the BLAST db
    if blast_files.len() < blast_suffixes.len() {
        info!("Making BLAST db");
        let mut makeblastdb_cmd = Command::new("makeblastdb");
        let makeblastdb_args = &[
            "-out".to_string(),
            args.consensi.to_string_lossy().to_string(),
            "-parse_seqids".to_string(),
            "-dbtype".to_string(),
            "nucl".to_string(),
            "-in".to_string(),
            args.consensi.to_string_lossy().to_string(),
        ];

        debug!("Running 'makeblastdb {}'", &makeblastdb_args.join(" "));
        let res = makeblastdb_cmd.args(makeblastdb_args).output()?;
        if !res.status.success() {
            bail!(String::from_utf8(res.stderr)?);
        }
    }

    let alignment_outfile = args.outdir.join("alignment.txt");
    let mut rmblastn_args = vec![
        "-num_alignments".to_string(),
        "9999999".to_string(),
        "-mask_level".to_string(),
        args.align_mask_level.to_string(),
        "-min_raw_gapped_score".to_string(),
        args.align_min_score.to_string(),
        "-num_threads".to_string(),
        args.num_threads.unwrap_or(num_cpus::get()).to_string(),
        "-db".to_string(),
        args.consensi.to_string_lossy().to_string(),
        "-query".to_string(),
        all_seqs_path.to_string_lossy().to_string(),
        "-out".to_string(),
        alignment_outfile.to_string_lossy().to_string(),
        "-outfmt".to_string(),
        "6 score qseqid sseqid".to_string(),
        // TODO: Verify these options
        "-gapopen".to_string(),
        "20".to_string(),
        "-gapextend".to_string(),
        "5".to_string(),
        "-word_size".to_string(),
        "7".to_string(),
        "-xdrop_ungap".to_string(),
        "400".to_string(),
        "-xdrop_gap_final".to_string(),
        "200".to_string(),
        "-xdrop_gap".to_string(),
        "100".to_string(),
        "-dust".to_string(),
        "no".to_string(),
        // TODO: verify remove
        //"-complexity_adjust".to_string(),
    ];

    let rmblastn = match &args.aligner {
        Some(path) => PathBuf::from(path.to_string()),
        _ => which("rmblastn")?,
    };
    let mut cmd = Command::new(&rmblastn);
    if let Some(perl5lib) = &args.perl5lib {
        cmd.env("PERL5LIB", perl5lib);
    }

    if let Some(matrix) = &args.alignment_matrix {
        let filename = matrix.file_name().unwrap_or_else(|| {
            panic!("Failed to get filename from '{}'", matrix.display())
        });

        let dirname = matrix.parent().unwrap_or_else(|| {
            panic!("Failed to get dirname from '{}'", matrix.display())
        });

        cmd.env("BLASTMAT", dirname);
        rmblastn_args.extend_from_slice(&[
            "-matrix".to_string(),
            filename.to_string_lossy().to_string(),
        ]);
    }

    debug!(
        "Running '{} {}'",
        rmblastn.display(),
        rmblastn_args.join(" ")
    );

    let start = Instant::now();
    let res = cmd.args(rmblastn_args).output()?;
    if !res.status.success() {
        bail!(String::from_utf8(res.stderr)?);
    }

    info!(
        "Alignment finished in {}",
        format_seconds(start.elapsed().as_secs())
    );

    Ok(alignment_outfile)
}

// --------------------------------------------------
fn call_winners(
    scores_file: &PathBuf,
    lambda: f64,
    confidence_margin: isize,
) -> Result<Winners> {
    debug!(
        r#"Calling winners from scores file "{}""#,
        scores_file.display()
    );

    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .from_path(scores_file)?;
    let records = reader.deserialize();

    let mut scores: HashMap<String, HashMap<String, u32>> = HashMap::new();

    for res in records {
        let rec: AlignmentScore = res?;
        scores
            .entry(rec.target)
            .or_default()
            .entry(rec.query)
            .and_modify(|v| *v = cmp::max(rec.score, *v))
            .or_insert(rec.score);
    }

    let mut clear_winners: HashMap<String, u32> = HashMap::new();
    let mut winning_sets: HashMap<StringPair, u32> = HashMap::new();
    for (_target, bit_scores) in scores.iter() {
        //debug!("\n>>> {target}");
        let families: Vec<_> = bit_scores.keys().collect();
        let raw_scores: Vec<_> = bit_scores.values().collect();
        let conf: Vec<_> = bitscore_to_confidence(&raw_scores, lambda)?;

        //if print {
        //    debug!("families   = {families:?}");
        //    debug!("raw scores = {raw_scores:?}");
        //    debug!(
        //        "confidence = {:?}",
        //        conf.iter()
        //            .map(|v| format!("{:0.04}", v))
        //            .collect::<Vec<_>>()
        //    );
        //}

        let pos: Vec<_> = (0..conf.len()).collect();
        let mut all_comps: Vec<bool> = vec![];
        for &i in &pos {
            let val = conf[i];
            let others: Vec<_> =
                pos.iter().filter(|&&j| j != i).map(|&j| conf[j]).collect();
            let pairs: Vec<_> = std::iter::repeat(val)
                .take(others.len())
                .zip(others)
                .collect();
            let comps: Vec<_> = pairs
                .iter()
                .map(|(x, y)| (x * (1. / confidence_margin as f64)) > *y)
                .collect();
            all_comps.push(comps.iter().all(|&v| v));
            //if print {
            //    debug!(
            //        "{:8} {:0.04} {:5} => {:?}",
            //        families[i],
            //        val,
            //        comps.iter().all(|&v| v),
            //        comps
            //    );
            //}
        }

        let winners: Vec<String> = families
            .iter()
            .zip(all_comps)
            .filter(|&(_, win)| win)
            .map(|(fam, _)| fam.to_string())
            .collect();

        if winners.len() == 1 {
            let winner = winners.first().unwrap();
            clear_winners
                .entry(winner.to_string())
                .and_modify(|v| *v += 1)
                .or_insert(1);
            //debug!("Clear Winner: {winner}");
        } else {
            let mut fam_comps: Vec<_> = conf.iter().zip(families).collect();
            fam_comps.sort_by(|a, b| {
                b.0.partial_cmp(a.0).unwrap().then_with(|| a.1.cmp(b.1))
            });
            let (&top_conf, _) = fam_comps.first().unwrap();
            let threshold = top_conf * (1. / confidence_margin as f64);
            let winning_set: Vec<_> = fam_comps
                .iter()
                .filter_map(|(&conf, fam)| (conf > threshold).then_some(fam))
                .collect();

            // The permutations will include A/B and B/A
            // It's important to store the symmetrical keys
            // Even though this is a duplication of the data
            for pair in winning_set.into_iter().permutations(2) {
                if let [&f1, &f2] = pair[..] {
                    let key = StringPair(f1.to_string(), f2.to_string());
                    winning_sets
                        .entry(key)
                        .and_modify(|v| *v += 1)
                        .or_insert(1);
                }
            }
        }
    }

    Ok(Winners {
        clear_winners,
        winning_sets,
    })
}

// --------------------------------------------------
fn cat_sequences(inputs: &[PathBuf], outpath: &PathBuf) -> Result<()> {
    let mut output = open_for_write(outpath)?;
    for filename in inputs.iter() {
        let mut reader = parse_reader(open(filename)?)?;
        let basename = Path::new(filename)
            .file_stem()
            .unwrap()
            .to_string_lossy()
            .to_string();

        while let Some(rec) = reader.iter_record()? {
            // TODO: Remove prefix of basename?
            writeln!(
                output,
                ">{}__{}{}\n{}",
                basename,
                rec.head(),
                rec.des(),
                rec.seq()
            )?;
        }
    }

    Ok(())
}

// --------------------------------------------------
fn check_family_instances(
    consensi_file: &PathBuf,
    instances_100: &[PathBuf],
) -> Result<HashMap<String, PathBuf>> {
    let mut reader = parse_reader(open(consensi_file)?)?;
    let mut consensi_names: HashMap<String, u32> = HashMap::new();
    while let Some(rec) = reader.iter_record()? {
        consensi_names
            .entry(rec.head().to_string())
            .and_modify(|v| *v += 1)
            .or_insert(1);
    }

    let mut dups: Vec<_> = consensi_names
        .iter()
        .flat_map(|(name, &count)| (count > 1).then_some(name))
        .collect();

    if !dups.is_empty() {
        // Have to sort for tests
        dups.sort();
        bail!(
            "The following consensi IDs are duplicated: {}",
            dups.iter().join(", ")
        );
    }

    let consensi_names: HashSet<String> =
        consensi_names.keys().cloned().collect();

    let mut family_to_path: HashMap<String, PathBuf> = HashMap::new();
    for instance in instances_100 {
        match instance.file_stem() {
            Some(stem) => family_to_path
                .insert(stem.to_string_lossy().to_string(), instance.clone()),
            _ => bail!(
                "Cannot get filename from {}",
                instance.to_string_lossy().to_string()
            ),
        };
    }

    let instance_names: HashSet<String> =
        family_to_path.keys().cloned().collect();

    let missing_instances: HashSet<_> =
        consensi_names.difference(&instance_names).collect();

    if !missing_instances.is_empty() {
        let mut missing_instances: Vec<_> =
            missing_instances.iter().collect();
        missing_instances.sort();
        bail!(
            "Missing the following instances: {}",
            missing_instances.iter().join(", ")
        );
    }

    let missing_consensi: HashSet<_> =
        instance_names.difference(&consensi_names).collect();

    if !missing_consensi.is_empty() {
        let mut missing_consensi: Vec<_> = missing_consensi.iter().collect();
        missing_consensi.sort();
        bail!(
            "Missing the following consensi: {}",
            missing_consensi.iter().join(", ")
        );
    }

    Ok(family_to_path)
}

// --------------------------------------------------
fn downsample(
    fasta: &PathBuf,
    num_wanted: usize,
    upper_range: usize,
    mut output: impl Write,
) -> Result<usize> {
    // Assumes FASTA file contains 100 longest sequences
    let mut reader = parse_reader(open(fasta)?)?;
    let mut rng = thread_rng();
    let range: Vec<usize> = (0..upper_range).collect();
    let mut take: Vec<usize> = range
        .choose_multiple(&mut rng, num_wanted)
        .cloned()
        .collect();
    take.sort();
    take.reverse();

    let mut next_take = take.pop().unwrap();
    let mut i = 0;
    let mut num_taken = 0;

    while let Some(rec) = reader.iter_record()? {
        if i == next_take {
            writeln!(output, ">{}{}\n{}", rec.head(), rec.des(), rec.seq())?;
            num_taken += 1;
            match take.pop() {
                // Keep skipping until we get to the next take value
                Some(val) => next_take = val,

                // We've run out of take values, so leave the loop
                _ => break,
            }
        }
        i += 1;
    }

    Ok(num_taken)
}

// --------------------------------------------------
fn extract_scores(
    alignment_file: &PathBuf,
    prev_scores_file: &Option<PathBuf>,
    consensi: &PathBuf,
    outdir: &Path,
) -> Result<PathBuf> {
    debug!(
        "Extracting scores from alignment '{}' {:?}",
        alignment_file.display(),
        prev_scores_file
    );

    let mut reader = parse_reader(open(consensi)?)?;
    let mut consensi_names: HashMap<String, String> = HashMap::new();
    while let Some(rec) = reader.iter_record()? {
        consensi_names.insert(
            rec.head().trim().to_string(),
            rec.des().trim().to_string(),
        );
    }

    let scores_file = outdir.join("alignment-scores.tsv");
    let mut scores_wtr = WriterBuilder::new()
        .has_headers(true)
        .delimiter(b'\t')
        .from_path(&scores_file)?;

    // This is a hash to remember the family names we handled
    // in this round. When reading the scores from the previous
    // round, we need to skip these queries.
    let mut skip_query: HashSet<String> = HashSet::new();

    //let file = open(alignment_file)?;
    // Regexes for parsing
    //let starts_with_score = Regex::new(r"^\d+\s").unwrap();
    //let spaces = Regex::new(r"\s+").unwrap();

    //for line in file.lines().map_while(Result::ok) {
    //    // Take only the score lines that start with an integer
    //    if starts_with_score.is_match(&line) {
    //        let parts: Vec<_> = spaces.split(&line).collect();
    //        if parts.len() == 12 {
    //            let score: u32 = parts[0].parse()?;
    //            let target = parts[4].to_string();
    //            let query = parts[8].to_string();

    //            match consensi_names.get(&query) {
    //                Some(query_name) => {
    //                    // Don't process this query again
    //                    for family in parse_newick(query_name) {
    //                        skip_query.insert(family);
    //                    }

    //                    scores_wtr.serialize(AlignmentScore {
    //                        score,
    //                        target,
    //                        query: query_name.clone(),
    //                        qseq: None,
    //                        sseq: None,
    //                    })?
    //                }
    //                _ => eprintln!("Cannot find query '{query}'"),
    //            }
    //        }
    //    }
    //}

    // BLAST output fails to include headers
    let mut align_reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(false)
        .from_path(alignment_file)?;
    let records = align_reader.records();
    for res in records {
        let res = res?;
        let rec: RmBlastOutput = res.deserialize(None)?;
        match consensi_names.get(rec.query) {
            Some(query_name) => {
                // Note the families involved in the previous round's
                // merges in order to skip them when adding the
                // previous round's scores.
                for family in parse_newick(query_name) {
                    skip_query.insert(family);
                }

                scores_wtr.serialize(AlignmentScore {
                    score: rec.score,
                    target: rec.target.to_string(),
                    query: query_name.clone(),
                })?
            }
            _ => eprintln!("Cannot find query '{}'", rec.query),
        }
    }

    // Add all the previous scores for queries not yet merged
    if let Some(prev_scores) = prev_scores_file {
        let mut reader = ReaderBuilder::new()
            .delimiter(b'\t')
            .has_headers(true)
            .from_path(prev_scores)?;
        let records = reader.deserialize();
        for res in records {
            let rec: AlignmentScore = res?;

            // Move along if any of this sequence's families have been seen
            let query_families = parse_newick(&rec.query);
            if !query_families.into_iter().any(|v| skip_query.contains(&v)) {
                scores_wtr.serialize(rec)?;
            }
        }
    }

    Ok(scores_file)
}

// --------------------------------------------------
fn find_independence(num_wins: u32, num_shared: u32) -> f64 {
    if num_shared > 0 {
        num_wins as f64 / (num_wins as f64 + num_shared as f64)
    } else {
        1.
    }
}

// --------------------------------------------------
/// Find the longest N number of sequences by grouping the sequences
/// by length and returning the shortest length that includes at least
/// the desired "number"
fn find_min_len(file: impl BufRead, number: usize) -> Result<usize> {
    let mut reader = parse_reader(file)?;
    let mut count_by_len: HashMap<usize, usize> = HashMap::new();
    while let Some(rec) = reader.iter_record()? {
        let len = rec.seq().len();
        count_by_len.entry(len).and_modify(|i| *i += 1).or_insert(1);
    }

    let mut lengths: Vec<&usize> = count_by_len.keys().collect();
    lengths.sort();
    lengths.reverse();

    let mut top_n = 0;
    let mut min_len = 0;
    for len in lengths {
        min_len = *len;
        top_n += count_by_len.get(len).unwrap_or(&0);
        if top_n >= number {
            break;
        }
    }

    Ok(min_len)
}

// --------------------------------------------------
fn format_seconds(seconds: u64) -> String {
    let mut delta = Duration::seconds(seconds as i64);
    let mut ret = vec![];
    let days = delta.num_days();
    if days > 0 {
        ret.push(format!("{days} day{}", if days == 1 { "" } else { "s" }));
    }
    delta -= Duration::seconds(days * 24 * 60 * 60);

    let hours = delta.num_hours();
    if hours > 0 {
        ret.push(format!(
            "{hours} hour{}",
            if hours == 1 { "" } else { "s" }
        ));
    }
    delta -= Duration::seconds(hours * 60 * 60);

    let minutes = delta.num_minutes();
    if minutes > 0 {
        ret.push(format!(
            "{minutes} minute{}",
            if minutes == 1 { "" } else { "s" }
        ));
    }
    delta -= Duration::seconds(minutes * 60);

    let seconds = delta.num_seconds();
    if seconds > 0 || ret.is_empty() {
        ret.push(format!(
            "{seconds} second{}",
            if seconds == 1 { "" } else { "s" }
        ));
    }

    ret.join(", ")
}

// --------------------------------------------------
fn independence(winners: Winners) -> Vec<Independence> {
    let mut families = HashSet::<&str>::new();
    for StringPair(f1, f2) in winners.winning_sets.keys() {
        families.insert(f1);
        families.insert(f2);
    }
    let mut vals = vec![];

    for pair in families.into_iter().permutations(2) {
        if let [f1, f2] = pair[..] {
            let num_wins = *winners.clear_winners.get(f1).unwrap_or(&0);
            let key = StringPair(f1.to_string(), f2.to_string());
            let &num_shared = winners.winning_sets.get(&key).unwrap_or(&0u32);
            let ind = find_independence(num_wins, num_shared);

            vals.push(Independence {
                f1: f1.to_string(),
                f2: f2.to_string(),
                val: ind,
            });
        }
    }

    // Sort from least independent to most
    vals.sort_by(|a, b| {
        a.val
            .partial_cmp(&b.val)
            .unwrap()
            .then_with(|| a.f1.cmp(&b.f1))
            .then_with(|| a.f2.cmp(&b.f2))
    });

    vals
}

// --------------------------------------------------
fn merge_families(
    args: MergeFamilies,
    family_to_path: &mut HashMap<String, PathBuf>,
) -> Result<String> {
    fs::create_dir_all(&args.outdir)?;

    let fams1 = parse_newick(&args.family1);
    let fams2 = parse_newick(&args.family2);
    let num_fams1 = fams1.len() as f64;
    let num_fams2 = fams2.len() as f64;
    let num_fams_total = num_fams1 + num_fams2;
    let num_seqs_total = 100;
    let num_from1 = (num_seqs_total as f64 * (num_fams1 / num_fams_total))
        .round() as usize;
    let num_from2 = (num_seqs_total as f64 * (num_fams2 / num_fams_total))
        .round() as usize;
    let f1 = fams1.join("::");
    let f2 = fams2.join("::");
    let new_family_name = format!("{f1}::{f2}");
    let new_family_path =
        args.instances_100_dir.join(format!("{new_family_name}.fa"));

    debug!(
        "Merging {num_from1} from {f1}, {num_from2} from {f2} => {}",
        new_family_path.display()
    );

    // Block to isolate "output" and force close when passing out of scope
    {
        let mut output = open_for_write(&new_family_path)?;
        let mut total_taken = 0;
        for (fam, num) in &[(f1, num_from1), (f2, num_from2)] {
            let fasta = family_to_path.get(fam).unwrap_or_else(|| {
                panic!("Missing instances for family '{fam}'")
            });
            total_taken +=
                downsample(fasta, *num, num_seqs_total, &mut output)?;
        }

        if total_taken == 0 {
            bail!("Failed to extract any sequences for MSA");
        }
    }

    // Remember this for the next iteration
    let _ = &family_to_path.insert(new_family_name, new_family_path.clone());

    // Copy the file to the current merge dir
    let msa_input = args.outdir.join("msa-input.fa");
    fs::copy(new_family_path, &msa_input)?;

    let refiner = match &args.refiner {
        Some(path) => PathBuf::from(path.to_string()),
        _ => which("Refiner")?,
    };

    let mut refiner_args =
        vec!["-threads".to_string(), args.num_threads.to_string()];

    if let Some(LogLevel::Debug) = args.log {
        refiner_args.push("-debug".to_string())
    }

    if let Some(rmblast_dir) = &args.rmblast_dir {
        refiner_args.extend_from_slice(&[
            "--rmblast_dir".to_string(),
            rmblast_dir.to_string(),
        ]);
    }

    refiner_args.push(msa_input.to_string_lossy().to_string());
    debug!(
        r#"Running "{} {}""#,
        &refiner.display(),
        &refiner_args.join(" ")
    );

    let start = Instant::now();
    let mut cmd = Command::new(&refiner);
    if let Some(perl5lib) = &args.perl5lib {
        cmd.env("PERL5LIB", perl5lib);
    }
    let res = cmd.args(refiner_args).output()?;
    if !res.status.success() {
        debug!("{}", String::from_utf8(res.stdout)?);
        bail!(String::from_utf8(res.stderr)?);
    }

    info!(
        "Refiner finished in {}",
        format_seconds(start.elapsed().as_secs())
    );

    let consensus_path = &args.outdir.join("msa-input.fa.refiner_cons");
    if !consensus_path.exists() {
        bail!(
            "Failed to find expected consensus file {}",
            consensus_path.display()
        );
    }

    let mut reader = parse_reader(open(consensus_path)?)?;
    let consensus_seq = reader
        .iter_record()?
        .map(|rec| rec.seq().to_string())
        .expect("Failed to read consensus file");

    Ok(consensus_seq)
}

// --------------------------------------------------
fn number_consensi(
    consensi: &PathBuf,
    outpath: &PathBuf,
) -> Result<HashMap<String, String>> {
    let mut outfile = open_for_write(outpath)?;
    let mut reader = parse_reader(open(consensi)?)?;
    let mut seqs: HashMap<String, String> = HashMap::new();
    let mut i = 0;

    while let Some(rec) = reader.iter_record()? {
        writeln!(outfile, ">{i} {}\n{}", rec.head(), rec.seq())?;
        seqs.insert(rec.head().to_string(), rec.seq().to_string());
        i += 1;
    }

    Ok(seqs)
}

// --------------------------------------------------
fn open(filename: &PathBuf) -> Result<Box<dyn BufRead>> {
    Ok(Box::new(BufReader::new(File::open(filename).map_err(
        |e| anyhow!("Cannot read {}: {e}", filename.display()),
    )?)))
}

// --------------------------------------------------
fn open_for_write(filename: &PathBuf) -> Result<Box<dyn Write>> {
    Ok(Box::new(BufWriter::new(File::create(filename).map_err(
        |e| anyhow!("Cannot write {}: {e}", filename.display()),
    )?)))
}

// --------------------------------------------------
fn parse_newick(val: &str) -> Vec<String> {
    let mut ret = vec![];
    match newick::from_string(val) {
        // Incoming value is only Newick when merged
        Ok(trees) => {
            let mut leaves: Vec<_> = trees
                .into_iter()
                .flat_map(|tree| {
                    let leaves = tree.leaves().collect::<Vec<_>>();
                    leaves.into_iter().flat_map(move |leaf| {
                        tree.name(leaf).map(|name| name.to_string())
                    })
                })
                .collect();
            ret.append(&mut leaves);
        }
        // Otherwise return the original string
        Err(_) => ret.push(val.to_string()),
    }

    ret
}

// --------------------------------------------------
fn take_longest_100(args: &Args, outdir: &Path) -> Result<Vec<PathBuf>> {
    // TODO: Make a parameter?
    let limit = 100;

    let mut ret = vec![];
    for file in &args.instances {
        let min_length = find_min_len(open(file)?, 100)?;
        let outfile = outdir.join(file.file_name().unwrap());
        let mut output = open_for_write(&outfile)?;
        let mut reader = parse_reader(open(file)?)?;
        let mut taken = 0;

        while let Some(rec) = reader.iter_record()? {
            if taken == limit {
                break;
            }

            let seq_len = rec.seq().len();
            if min_length > 0 && seq_len >= min_length {
                taken += 1;
                writeln!(
                    output,
                    ">{}{}\n{}",
                    rec.head(),
                    rec.des(),
                    rec.seq()
                )?;
            }
        }

        ret.push(outfile);
    }

    Ok(ret)
}

// --------------------------------------------------
#[cfg(test)]
mod tests {
    use super::{
        align, bitscore_to_confidence, call_winners, cat_sequences,
        check_family_instances, downsample, extract_scores,
        find_independence, find_min_len, format_seconds, independence,
        number_consensi, open, parse_newick, Args, Independence, StringPair,
        Winners,
    };
    use anyhow::Result;
    use kseq::parse_reader;
    use pretty_assertions::assert_eq;
    use std::{
        collections::HashMap,
        fs::{self, File},
        io::Cursor,
        path::PathBuf,
    };
    use tempfile::{tempdir, NamedTempFile};

    // Or "tests/inputs/25p41g.matrix"?
    const MATRIX: &str =
        "/Users/kyclark/work/RepeatMasker/Matrices/ncbi/nt/25p41g.matrix";

    #[test]
    fn test_align() -> Result<()> {
        let outdir = tempdir()?;
        let args = Args {
            consensi: PathBuf::from("tests/inputs/consensi.fa".to_string()),
            instances: vec![PathBuf::from(
                "tests/inputs/test_set.fa".to_string(),
            )],
            alignment_matrix: Some(MATRIX.into()),
            outdir: outdir.path().into(),
            perl5lib: Some("/Users/kyclark/work/RepeatMasker".to_string()),
            lambda: 0.1227,
            confidence_margin: 3,
            independence_threshold: 0.5,
            aligner: Some(
                "/Users/kyclark/work/RepeatModeler/util/align.pl".to_string(),
            ),
            refiner: Some(
                "/Users/kyclark/work/RepeatModeler/Refiner".to_string(),
            ),
            rmblast_dir: Some("/Users/kyclark/.local/bin".to_string()),
            num_threads: None,
            align_gap_init: -25,
            align_extension: -5,
            align_min_match: 7,
            align_bandwidth: 14,
            align_mask_level: 101,
            align_min_score: 200,
            log: None,
            log_file: None,
        };
        let all_seqs_path = PathBuf::from("tests/inputs/all_seqs.fa");
        let res = align(&args, &all_seqs_path);
        assert!(res.is_ok());

        let alignment_file = res.unwrap();
        assert!(alignment_file.exists());

        // Output seems non-deterministic?
        //let actual = fs::read_to_string(alignment_file)?;
        //let expected = fs::read_to_string("tests/outputs/alignment.txt")?;
        //assert_eq!(actual, expected);
        Ok(())
    }

    #[test]
    fn test_bitscore_to_confidence() -> Result<()> {
        let res = bitscore_to_confidence(
            &[&2274, &2234, &2224, &2245, &2295],
            0.1227,
        );
        assert!(res.is_ok());
        assert_eq!(
            res.unwrap(),
            [
                0.14088157852127628,
                0.004692442832468759,
                0.002004634432314714,
                0.011959118673026447,
                0.8404622255409138,
            ]
        );

        Ok(())
    }

    #[test]
    fn test_call_winners() -> Result<()> {
        let scores_file = PathBuf::from("tests/outputs/alignment-scores.tsv");
        let lambda = 0.1227;
        let confidence_margin = 3;
        let res = call_winners(&scores_file, lambda, confidence_margin);
        assert!(res.is_ok());

        let winners = res.unwrap();
        let expected_wins: Vec<(&str, u32)> = vec![
            ("AluY", 95),
            ("AluYm1", 72),
            ("AluYb9", 93),
            ("AluYb8", 103),
            ("AluYa5", 98),
        ];

        for (key, val) in expected_wins.iter() {
            assert_eq!(winners.clear_winners.get(*key), Some(val));
        }

        // Winning sets should be symmetrical (A/B, B/A)
        let expected_winning_sets: HashMap<StringPair, u32> =
            HashMap::from([
                (StringPair("AluYm1".to_string(), "AluY".to_string()), 28),
                (StringPair("AluY".to_string(), "AluYm1".to_string()), 28),
                //
                (StringPair("AluYb8".to_string(), "AluYa5".to_string()), 3),
                (StringPair("AluYa5".to_string(), "AluYb8".to_string()), 3),
                //
                (StringPair("AluYb8".to_string(), "AluYb9".to_string()), 4),
                (StringPair("AluYb9".to_string(), "AluYb8".to_string()), 4),
                //
                (StringPair("AluY".to_string(), "AluYb8".to_string()), 3),
                (StringPair("AluYb8".to_string(), "AluY".to_string()), 3),
                //
                (StringPair("AluY".to_string(), "AluYa5".to_string()), 7),
                (StringPair("AluYa5".to_string(), "AluY".to_string()), 7),
            ]);

        for (pair, score) in expected_winning_sets {
            let res = winners.winning_sets.get(&pair);
            assert!(res.is_some());
            assert_eq!(res.unwrap(), &score);
        }

        Ok(())
    }

    #[test]
    fn test_cat_sequences() -> Result<()> {
        let inputs = vec![
            PathBuf::from("tests/inputs/instances_100/AluY.fa"),
            PathBuf::from("tests/inputs/instances_100/AluYa5.fa"),
            PathBuf::from("tests/inputs/instances_100/AluYb8.fa"),
            PathBuf::from("tests/inputs/instances_100/AluYb9.fa"),
            PathBuf::from("tests/inputs/instances_100/AluYm1.fa"),
        ];
        let outdir = tempdir()?;
        let outpath = outdir.path().join("all_seqs.fa");
        let res = cat_sequences(&inputs, &outpath);
        assert!(res.is_ok());
        assert!(outpath.exists());

        let mut reader = parse_reader(open(&outpath)?)?;
        let mut count = 0;
        while (reader.iter_record()?).is_some() {
            count += 1;
        }
        assert_eq!(count, 500);

        Ok(())
    }

    #[test]
    fn test_check_family_instances() -> Result<()> {
        // The consensi contains duplicated IDs
        let consensi = PathBuf::from("tests/inputs/dup_consensi.fa");
        let instances_100 = &[PathBuf::from(
            "tests/inputs/instances_100/AluY.fa".to_string(),
        )];
        let res = check_family_instances(&consensi, instances_100);
        assert!(res.is_err());
        assert_eq!(
            res.unwrap_err().to_string(),
            "The following consensi IDs are duplicated: AluYb9, AluYm1"
        );

        // The consensi has 5 unique IDs but the instances only 2
        let consensi = PathBuf::from("tests/inputs/consensi.fa");
        let instances_100 = &[
            PathBuf::from("tests/inputs/AluY.fa".to_string()),
            PathBuf::from("tests/inputs/AluYm1.fa".to_string()),
        ];
        let res = check_family_instances(&consensi, instances_100);
        assert!(res.is_err());
        assert_eq!(
            res.unwrap_err().to_string(),
            "Missing the following instances: AluYa5, AluYb8, AluYb9",
        );

        // The consensi has 2 unique IDs but the instances have 4
        let consensi = PathBuf::from("tests/inputs/two_consensi.fa");
        let instances_100 = vec![
            PathBuf::from("tests/inputs/AluY.fa".to_string()),
            PathBuf::from("tests/inputs/AluYa5.fa".to_string()),
            PathBuf::from("tests/inputs/AluYb9.fa".to_string()),
            PathBuf::from("tests/inputs/AluYm1.fa".to_string()),
        ];
        let res = check_family_instances(&consensi, &instances_100);
        assert!(res.is_err());
        assert_eq!(
            res.unwrap_err().to_string(),
            "Missing the following consensi: AluYb9, AluYm1",
        );

        Ok(())
    }

    #[test]
    fn test_downsample() -> Result<()> {
        let fasta = PathBuf::from("tests/inputs/AluY_5.fa");
        let tmp_dir = tempdir()?;
        let outpath = tmp_dir.path().join("sub.fa");

        // Scoped to force close of output filehandle
        {
            let out = File::create(&outpath)?;

            // Select 3 of the 5 sequences
            let res = downsample(&fasta, 3, 5, &out);
            assert!(res.is_ok());
            assert_eq!(res.unwrap(), 3);
        }

        {
            let sub = File::open(&outpath)?;
            let mut reader = parse_reader(sub)?;
            let mut num = 0;
            while (reader.iter_record()?).is_some() {
                num += 1;
            }

            // Ensure 3 sequences were written
            assert_eq!(num, 3);
        }

        // Try to select more than the 5 sequences
        // Should only get the actual 5
        {
            let out = File::create(&outpath)?;
            let res = downsample(&fasta, 10, 5, &out);
            assert!(res.is_ok());
            assert_eq!(res.unwrap(), 5);
        }

        {
            let sub = File::open(&outpath)?;
            let mut reader = parse_reader(sub)?;
            let mut num = 0;
            while (reader.iter_record()?).is_some() {
                num += 1;
            }

            // Ensure 5 sequences were written
            assert_eq!(num, 5);
        }

        Ok(())
    }

    #[test]
    fn test_extract_scores() -> Result<()> {
        let outdir = tempdir()?;
        let consensi =
            PathBuf::from("tests/inputs/numbered_consensi.fa".to_string());
        let alignment_file = PathBuf::from("tests/inputs/alignment.txt");
        let orig_scores: Option<PathBuf> = None;
        let res = extract_scores(
            &alignment_file,
            &orig_scores,
            &consensi,
            outdir.path(),
        );
        assert!(res.is_ok());

        let scores_file = res.unwrap();
        assert!(scores_file.exists());

        let actual = fs::read_to_string(scores_file)?;
        let expected =
            fs::read_to_string("tests/outputs/alignment-scores.tsv")?;
        assert_eq!(actual, expected);
        Ok(())
    }

    #[test]
    fn test_find_independence() -> Result<()> {
        // No shared wins means fully independent
        //                       wins  shared
        let res = find_independence(1, 0);
        assert_eq!(res, 1.);

        // One clear winner, one shared win == 50%
        //                       wins  shared
        let res = find_independence(1, 1);
        assert_eq!(res, 0.5);

        // No clear winner, only shared wins == 0%
        //                       wins  shared
        let res = find_independence(0, 2);
        assert_eq!(res, 0.);

        //                        wins  shared
        let res = find_independence(37, 46);
        assert_eq!(res, 0.4457831325301205);
        Ok(())
    }

    #[test]
    fn test_find_min_len() -> Result<()> {
        let fasta = ">1\nAC\n>2\nA\n>3\nACG\n>4\nACGT\n>5\nACGTAA\n>6\nAC\n";

        // Asking for more sequences (10) than present (6) returns the
        // shortest sequence length (1)
        let res = find_min_len(Cursor::new(fasta), 10);
        assert!(res.is_ok());
        assert_eq!(res.unwrap(), 1);

        // Asking for one sequence finds the longest (6bp)
        let res = find_min_len(Cursor::new(fasta), 1);
        assert!(res.is_ok());
        assert_eq!(res.unwrap(), 6);

        let res = find_min_len(Cursor::new(fasta), 4);
        assert!(res.is_ok());
        assert_eq!(res.unwrap(), 2);
        Ok(())
    }

    #[test]
    fn test_format_seconds() -> Result<()> {
        let one_hour = 60 * 60;
        let one_day = one_hour * 24;
        assert_eq!(format_seconds(0), "0 seconds");
        assert_eq!(format_seconds(1), "1 second");
        assert_eq!(format_seconds(59), "59 seconds");
        assert_eq!(format_seconds(60), "1 minute");
        assert_eq!(format_seconds(120), "2 minutes");
        assert_eq!(format_seconds(121), "2 minutes, 1 second");
        assert_eq!(format_seconds(one_hour), "1 hour");
        assert_eq!(format_seconds(one_hour + 1), "1 hour, 1 second");
        assert_eq!(
            format_seconds(one_hour + 121),
            "1 hour, 2 minutes, 1 second"
        );
        assert_eq!(
            format_seconds((one_hour * 4) + 59),
            "4 hours, 59 seconds"
        );
        assert_eq!(format_seconds(one_day), "1 day");
        assert_eq!(format_seconds(one_day + 2), "1 day, 2 seconds");
        Ok(())
    }

    #[test]
    fn test_independence() -> Result<()> {
        let clear_winners = HashMap::from([
            ("AluY".to_string(), 37),
            ("AluYm1".to_string(), 36),
            ("AluYb8".to_string(), 50),
        ]);

        // Winning sets are expected to have symmetrical keys
        let winning_sets = HashMap::from([
            (StringPair("AluYa5".to_string(), "AluYb8".to_string()), 6),
            (StringPair("AluYb8".to_string(), "AluYa5".to_string()), 6),
            //
            (StringPair("AluYa5".to_string(), "AluYm1".to_string()), 4),
            (StringPair("AluYm1".to_string(), "AluYa5".to_string()), 4),
            //
            (StringPair("AluY".to_string(), "AluYb8".to_string()), 6),
            (StringPair("AluYb8".to_string(), "AluY".to_string()), 6),
            //
            (StringPair("AluY".to_string(), "AluYa5".to_string()), 12),
            (StringPair("AluYa5".to_string(), "AluY".to_string()), 12),
            //
            (StringPair("AluY".to_string(), "AluYm1".to_string()), 46),
            (StringPair("AluYm1".to_string(), "AluY".to_string()), 46),
        ]);

        let ind = independence(Winners {
            clear_winners,
            winning_sets,
        });

        let expected = [
            Independence {
                f1: "AluYa5".to_string(),
                f2: "AluY".to_string(),
                val: 0.0,
            },
            Independence {
                f1: "AluYa5".to_string(),
                f2: "AluYb8".to_string(),
                val: 0.0,
            },
            Independence {
                f1: "AluYa5".to_string(),
                f2: "AluYm1".to_string(),
                val: 0.0,
            },
            Independence {
                f1: "AluYm1".to_string(),
                f2: "AluY".to_string(),
                val: 0.43902439024390244,
            },
            Independence {
                f1: "AluY".to_string(),
                f2: "AluYm1".to_string(),
                val: 0.4457831325301205,
            },
            Independence {
                f1: "AluY".to_string(),
                f2: "AluYa5".to_string(),
                val: 0.7551020408163265,
            },
            Independence {
                f1: "AluY".to_string(),
                f2: "AluYb8".to_string(),
                val: 0.8604651162790697,
            },
            Independence {
                f1: "AluYb8".to_string(),
                f2: "AluY".to_string(),
                val: 0.8928571428571429,
            },
            Independence {
                f1: "AluYb8".to_string(),
                f2: "AluYa5".to_string(),
                val: 0.8928571428571429,
            },
            Independence {
                f1: "AluYm1".to_string(),
                f2: "AluYa5".to_string(),
                val: 0.9,
            },
            Independence {
                f1: "AluYb8".to_string(),
                f2: "AluYm1".to_string(),
                val: 1.0,
            },
            Independence {
                f1: "AluYm1".to_string(),
                f2: "AluYb8".to_string(),
                val: 1.0,
            },
        ];
        assert_eq!(ind, expected);
        Ok(())
    }

    #[test]
    fn test_number_consensi() -> Result<()> {
        let orig_consensi = PathBuf::from("tests/inputs/consensi.fa");
        let outpath = NamedTempFile::new()?;
        let res = number_consensi(&orig_consensi, &outpath.path().into());
        assert!(res.is_ok());

        let actual = fs::read_to_string(outpath)?;
        let expected =
            fs::read_to_string("tests/outputs/numbered_consensi.fa")?;
        assert_eq!(actual, expected);

        Ok(())
    }

    #[test]
    fn test_parse_newick() -> Result<()> {
        let res = parse_newick("A");
        assert_eq!(res, vec!["A"]);

        let res = parse_newick("(((B,D):0.1,A):0.3,C):0.2;");
        assert_eq!(res, vec!["B", "D", "A", "C"]);

        let res = parse_newick("((X,Y):0.04,Z):0.25;");
        assert_eq!(res, vec!["X", "Y", "Z"]);
        Ok(())
    }
}
