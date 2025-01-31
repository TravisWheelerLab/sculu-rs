mod graph;

use anyhow::{anyhow, bail, Result};
use assert_cmd::Command;
use chrono::Duration;
use clap::Parser;
use csv::{ReaderBuilder, WriterBuilder};
use itertools::Itertools;
use kseq::parse_reader;
use log::debug;
use newick::Newick;
use noodles_fasta::{
    self,
    io::Reader as FastaReader,
    io::Writer as FastaWriter,
    record::{Definition as FastaDefinition, Record as FastaRecord},
};
use rand::{seq::SliceRandom, thread_rng};
use regex::Regex;
use serde::{Deserialize, Serialize};
use std::{
    cmp,
    collections::{HashMap, HashSet},
    fmt,
    fs::{self, File},
    io::{BufRead, BufReader, BufWriter, Write},
    path::{Path, PathBuf},
    time::Instant,
};
use tempfile::tempdir;
//use tabled::{Table, Tabled};
//use walkdir::WalkDir;
use which::which;

/// SCULU subfamily clustering tool
#[derive(Debug, Parser, Clone)]
#[command(author, version, about)]
pub struct Args {
    /// FASTA file of subfamily consensi
    #[arg(short, long, value_name = "CONSENSI")]
    pub consensi: PathBuf,

    /// Instance files for each subfamily
    #[arg(short, long, value_name = "INSTANCES", required = true, num_args=1..)]
    pub instances: Vec<PathBuf>,

    /// Output file
    #[arg(short, long, value_name = "OUTFILE", default_value = "families.fa")]
    pub outfile: PathBuf,

    /// Output directory (if you want working files preserved)
    #[arg(long, value_name = "OUTDIR")]
    pub outdir: Option<PathBuf>,

    /// Lambda value
    #[arg(long, value_name = "LAMBDA", default_value = "0.1227")]
    pub lambda: f64,

    /// Independence threshold
    #[arg(long, value_name = "IND", default_value = "0.5")]
    pub independence_threshold: f64,

    /// Confidence margin
    #[arg(long, value_name = "CONF", default_value = "3")]
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

    /// Alignment matrix
    #[arg(long, value_name = "MATRIX")]
    pub align_matrix: Option<PathBuf>,

    /// Alignment gap open penalty
    #[arg(long, value_name = "ALIGN_GAP_OPEN", default_value = "20")]
    pub align_gap_open: usize,

    /// Alignment gap extension penalty
    #[arg(long, value_name = "ALIGN_GAP_EXT", default_value = "5")]
    pub align_gap_extension: usize,

    /// Alignment word size
    #[arg(long, value_name = "ALIGN_WORD_SIZE", default_value = "7")]
    pub align_word_size: usize,

    /// Alignment x-drop gap
    #[arg(long, value_name = "ALIGN_XDROP_GAP", default_value = "100")]
    pub align_xdrop_gap: usize,

    /// Alignment x-drop ungap
    #[arg(long, value_name = "ALIGN_XDROP_UNGAP", default_value = "400")]
    pub align_xdrop_ungap: usize,

    /// Alignment x-drop final
    #[arg(long, value_name = "ALIGN_XDROP_FINAL", default_value = "200")]
    pub align_xdrop_final: usize,

    /// Alignment dust option
    #[arg(long)]
    pub align_dust: bool,

    /// Alignment complexity adjust
    #[arg(long)]
    pub align_complexity_adjust: bool,

    /// Alignment mask level
    #[arg(long, value_name = "MASKLEVEL", default_value = "101")]
    pub align_mask_level: i64,

    /// Alignment minimum score
    #[arg(long, value_name = "MINSCORE", default_value = "200")]
    pub align_min_score: i64,
}

#[derive(Debug, Deserialize, Clone)]
struct RmBlastOutput {
    score: u32,
    target: String,
    query: String,
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
    instances_dir: &'a PathBuf,
    num_threads: usize,
    refiner: &'a Option<String>,
    perl5lib: &'a Option<String>,
}

// --------------------------------------------------
pub fn run(args: Args) -> Result<()> {
    if let Some(path) = &args.align_matrix {
        if !path.is_file() {
            bail!("--align-matrix '{}' does not exist", path.display());
        }
    }

    // Final output file
    if args.outfile.is_file() && args.outfile.metadata()?.len() > 0 {
        bail!(r#"--outfile "{}" exists"#, args.outfile.display())
    }
    let outfile = BufWriter::new(File::create(&args.outfile)?);
    let mut fasta_writer = FastaWriter::new(outfile);

    // Use a tempdir if no named outdir
    let outdir = args
        .outdir
        .clone()
        .unwrap_or(tempdir()?.path().to_path_buf());
    if !outdir.is_dir() {
        fs::create_dir_all(&outdir)?;
    }

    // All logging goes into outdir
    env_logger::Builder::new()
        .filter_level(log::LevelFilter::Debug)
        .target(env_logger::Target::Pipe(Box::new(BufWriter::new(
            File::create(outdir.join("debug.log"))?,
        ))))
        .init();
    debug!("args = {:#?}", &args);

    // Create a directory to hold sampled instances
    // Copy the given instances to the working directory
    //
    // At first, I used all the instances directly
    // Then we decided to try downsampling methods like
    // taking the longest 100 or those longer than 50%
    // of the consensi, but now we're back to using all
    // and it's kind of silly to copy them like this.
    //
    // TODO: Go back to using the original Vec<PathBuf>?
    let instances_dir = &outdir.join("instances");
    debug!("Writing instances to '{}'", instances_dir.display());
    fs::create_dir_all(instances_dir)?;
    for file in &args.instances {
        fs::copy(
            file,
            instances_dir.join(file.file_name().expect("filename")),
        )?;
    }

    // Align consensi to themselves to find groups/batches
    let mut components = align_consensi_to_self(&outdir, &args)?;
    components.sort_by_key(|v| v.len());
    debug!("After aligning consensi to self, components =\n{components:#?}");

    let mut batch_num = 0;
    for families in components {
        if families.len() == 1 {
            debug!("Extracting independent family '{}'", families[0]);
            copy_fasta(&families, &args.consensi, &mut fasta_writer)?
        } else {
            batch_num += 1;
            run_batch(
                &outdir,
                batch_num,
                families,
                args.clone(),
                instances_dir,
                &mut fasta_writer,
            )?
        }
    }

    println!(
        r#"Done, see final families in "{}""#,
        args.outfile.display()
    );
    Ok(())
}

// --------------------------------------------------
fn copy_fasta<W: Write>(
    family_names: &[String],
    source: &PathBuf,
    destination: &mut FastaWriter<W>,
) -> Result<()> {
    let mut reader = FastaReader::new(BufReader::new(File::open(source)?));

    for result in reader.records() {
        let record = result?;
        let record_name = String::from_utf8(record.name().to_vec())?;
        if family_names.contains(&record_name) {
            destination.write_record(&record)?;
        }
    }

    Ok(())
}

// --------------------------------------------------
fn run_batch<W: Write>(
    outdir: &Path,
    batch_num: usize,
    families: Vec<String>,
    args: Args,
    instances_dir: &PathBuf,
    final_writer: &mut FastaWriter<W>,
) -> Result<()> {
    debug!(
        "==> Running batch {batch_num} for families {} <==",
        families.join(", ")
    );
    let start = Instant::now();
    let batch_dir = outdir.join(format!("batch{batch_num:03}"));
    fs::create_dir_all(&batch_dir)?;

    // Create a consensi file for this round containing only the given families
    let mut consensi = batch_dir.join("consensi.fa");
    {
        // Scoped to cause fasta_writer to close
        let outfile = BufWriter::new(File::create(&consensi)?);
        let mut fasta_writer = FastaWriter::new(outfile);
        copy_fasta(&families, &args.consensi, &mut fasta_writer)?;
    }

    // Check that the "consensi" sequences have clear relationships
    // to "instance" files. Returns a mutable hashmap from the
    // family names to their instance files for merges.
    debug!("Checking family names/instances");
    let mut family_to_path = check_family_instances(&consensi, instances_dir)?;
    debug!("family_to_path = '{family_to_path:#?}'");

    // Create query file for families
    let all_seqs_path = &batch_dir.join("all_seqs.fa");
    debug!(
        "Concatenating all instance sequences into '{}'",
        all_seqs_path.display()
    );

    // Put the instances for the given families into a single file
    cat_sequences(&args.instances, &families, all_seqs_path)?;

    // Create a starting directory for the merge iterations.
    let round0_dir = batch_dir.join("round000");
    fs::create_dir_all(&round0_dir)?;

    // As the consensi are merged, the family names are concatenated
    // into Newick-formatted strings to track the merges.
    // These can get quite long, causing `makeblastdb` to fail.
    // Make the sequence IDs simple integers by numbering
    // and move the family names to the description.
    let consensi_path = round0_dir.join("consensi.fa");
    debug!("Writing numbered consensi to '{}'", consensi_path.display());
    let mut consensi_seqs = number_consensi(&consensi, &consensi_path)?;
    consensi = consensi_path;

    // Start merging
    let trailing_semi = Regex::new(";$").unwrap();
    let mut prev_scores: Option<PathBuf> = None;
    for round in 1.. {
        debug!(">>> Round {round} <<<");
        let round_dir = batch_dir.join(format!("round{round:03}"));
        fs::create_dir_all(&round_dir)?;

        // Align all the sequences to the current consensi.
        // On the first round, all the original consensi will be included.
        // On future rounds, only the newly merged consensi will be present.
        let alignment_file = run_rmblastn(&args, &round_dir, &consensi, all_seqs_path)?;

        // Extract the scores from the alignment file.
        // On the first round, there will be no "prev_scores" file.
        // On the following rounds, the previous round's scores will be
        // included for those families that were not merged.
        let scores_file =
            extract_scores(&alignment_file, &prev_scores, &consensi, &round_dir)?;

        // Save this round's alignment scores for the next
        prev_scores = Some(scores_file.clone());

        // Find the clear/ambiguous winners from the scores
        let winners = call_winners(&scores_file, args.lambda, args.confidence_margin)?;

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
        let new_consensi_path = &round_dir.join("new-consensi.fa");
        let mut new_consensi = open_for_write(new_consensi_path)?;
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

            debug!(
                "{}: Merge {} => {} Ind: {:0.04}{other_score}",
                merge_num, &pair.f1, &pair.f2, pair.val
            );

            // Place all merge artefacts into a directory.
            let msa_dir = &round_dir.join(format!("msa-{merge_num:02}"));
            let new_consensus_seq = merge_families(
                MergeFamilies {
                    family1: pair.f1.clone(),
                    family2: pair.f2.clone(),
                    outdir: msa_dir.to_path_buf(),
                    instances_dir,
                    num_threads: args.num_threads.unwrap_or(num_cpus::get()),
                    refiner: &args.refiner,
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

            let f1_len = consensi_seqs.get(&pair.f1).map_or(0, |v| v.len());
            let f2_len = consensi_seqs.get(&pair.f2).map_or(0, |v| v.len());
            let new_len = new_consensus_seq.len();

            debug!(
                "{} len was {f1_len}, {} len was {f2_len}, \
                    new consensi len is {new_len}",
                pair.f1, pair.f2
            );

            // Update the consensi mapping
            consensi_seqs.remove(&pair.f1);
            consensi_seqs.remove(&pair.f2);
            consensi_seqs.insert(new_family_newick, new_consensus_seq.clone());

            // Keep track of all the families that have been merged
            for family in all_families {
                already_merged.insert(family);
            }
        }

        // Write the merged families to the new consensi file
        debug!(
            "Writing {} consensi to '{}'",
            &consensi_seqs.len(),
            new_consensi_path.display()
        );

        for (family_number, (family_name, seq)) in consensi_seqs.iter().enumerate() {
            writeln!(new_consensi, ">{family_number} {family_name}\n{seq}",)?;
        }

        consensi = new_consensi_path.to_path_buf();
    }

    let mut reader = FastaReader::new(BufReader::new(File::open(consensi)?));
    let mut new_seqs = 0;
    for result in reader.records() {
        let mut record = result?;
        if let Some(desc) = record.description() {
            record = FastaRecord::new(
                FastaDefinition::new(desc, None),
                record.sequence().clone(),
            )
        }
        final_writer.write_record(&record)?;
        new_seqs += 1;
    }

    debug!(
        "Added {new_seqs} famil{} in {}.",
        if new_seqs == 1 { "y" } else { "ies" },
        format_seconds(start.elapsed().as_secs()),
    );

    Ok(())
}

// --------------------------------------------------
// Align the consensi to itself to identify clusters
// [
//     [ "Charlie13a" ],
//     [ "Tigger3c" ],
//     [ "AluSc", "AluYm1", "AluYh3", "AluYh9", "AluYb8", "AluYb9" ],
//     [ "Charlie1a", "Charlie1", "Charlie2a", "Charlie2b" ],
// ]
fn align_consensi_to_self(outdir: &Path, args: &Args) -> Result<Vec<Vec<String>>> {
    let outdir = outdir.join("consensi_cluster");
    let output = run_rmblastn(args, &outdir, &args.consensi, &args.consensi)?;
    Ok(graph::connected_components(parse_alignment(&output)?))
}

// --------------------------------------------------
fn parse_alignment(blast_out: &PathBuf) -> Result<Vec<RmBlastOutput>> {
    // BLAST output fails to include headers
    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(false)
        .from_path(blast_out)?;

    let mut records = vec![];
    for res in reader.records() {
        let rec = res?;
        let blast_rec: RmBlastOutput = rec.deserialize(None)?;
        records.push(blast_rec);
    }

    Ok(records)
}

// --------------------------------------------------
fn run_rmblastn(
    args: &Args,
    outdir: &PathBuf,
    db: &Path,
    query: &Path,
) -> Result<PathBuf> {
    fs::create_dir_all(outdir)?;
    let db_path = &outdir.join("db");
    let makeblastdb_args = &[
        "-out".to_string(),
        db_path.to_string_lossy().to_string(),
        "-parse_seqids".to_string(),
        "-dbtype".to_string(),
        "nucl".to_string(),
        "-in".to_string(),
        db.to_string_lossy().to_string(),
    ];

    debug!("Running 'makeblastdb {}'", &makeblastdb_args.join(" "));
    let makeblastdb = which("makeblastdb").map_err(|e| anyhow!("makeblastdb: {e}"))?;
    let res = Command::new(makeblastdb).args(makeblastdb_args).output()?;

    if !res.status.success() {
        bail!(String::from_utf8(res.stderr)?);
    }

    let outfile = outdir.join("blast.out");
    let mut rmblastn_args = vec![
        "-db".to_string(),
        db_path.to_string_lossy().to_string(),
        "-query".to_string(),
        query.to_string_lossy().to_string(),
        "-out".to_string(),
        outfile.to_string_lossy().to_string(),
        "-outfmt".to_string(),
        "6 score qseqid sseqid".to_string(),
        "-num_threads".to_string(),
        args.num_threads.unwrap_or(num_cpus::get()).to_string(),
        "-num_alignments".to_string(),
        "9999999".to_string(),
        "-mask_level".to_string(),
        args.align_mask_level.to_string(),
        "-min_raw_gapped_score".to_string(),
        args.align_min_score.to_string(),
        "-gapopen".to_string(),
        args.align_gap_open.to_string(),
        "-gapextend".to_string(),
        args.align_gap_extension.to_string(),
        "-word_size".to_string(),
        args.align_word_size.to_string(),
        "-xdrop_ungap".to_string(),
        args.align_xdrop_ungap.to_string(),
        "-xdrop_gap_final".to_string(),
        args.align_xdrop_final.to_string(),
        "-xdrop_gap".to_string(),
        args.align_xdrop_gap.to_string(),
        "-dust".to_string(),
        if args.align_dust {
            "yes".to_string()
        } else {
            "no".to_string()
        },
    ];

    if args.align_complexity_adjust {
        rmblastn_args.push("-complexity_adjust".to_string());
    }

    let rmblastn = match &args.aligner {
        Some(path) => PathBuf::from(path.to_string()),
        _ => which("rmblastn").map_err(|e| anyhow!("rmblastn: {e}"))?,
    };

    let mut cmd = Command::new(&rmblastn);
    if let Some(perl5lib) = &args.perl5lib {
        cmd.env("PERL5LIB", perl5lib);
    }

    if let Some(matrix) = &args.align_matrix {
        let matrix_filename = matrix.file_name().unwrap_or_else(|| {
            panic!("Failed to get filename from '{}'", matrix.display())
        });

        let matrix_dir = matrix.parent().unwrap_or_else(|| {
            panic!("Failed to get dirname from '{}'", matrix.display())
        });

        cmd.env("BLASTMAT", matrix_dir);
        rmblastn_args.extend_from_slice(&[
            "-matrix".to_string(),
            matrix_filename.to_string_lossy().to_string(),
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

    debug!(
        "Rmblastn finished in {}",
        format_seconds(start.elapsed().as_secs())
    );

    Ok(outfile)
}

// --------------------------------------------------
fn bitscore_to_confidence(vals: &[&u32], lambda: f64) -> Result<Vec<f64>> {
    let mut converted: Vec<_> =
        vals.iter().map(|&&v| (v as f64 * lambda).exp2()).collect();

    if converted.iter().any(|v| v.is_infinite()) {
        // Scale all numbers down
        let delta = **vals.iter().max().unwrap() as i64 - 500;
        converted = vals
            .iter()
            .map(|&&v| ((v as i64 - delta) as f64 * lambda).exp2())
            .collect::<Vec<_>>();
    }

    let total: f64 = converted.iter().sum();
    if total > 0. {
        Ok(converted.into_iter().map(|v| v / total).collect())
    } else {
        bail!("Sum of converted values equals zero")
    }
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
        let families: Vec<_> = bit_scores.keys().collect();
        let raw_scores: Vec<_> = bit_scores.values().collect();
        let conf: Vec<_> = bitscore_to_confidence(&raw_scores, lambda)?;
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
                    winning_sets.entry(key).and_modify(|v| *v += 1).or_insert(1);
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
fn cat_sequences(
    inputs: &[PathBuf],
    families: &[String],
    outpath: &PathBuf,
) -> Result<()> {
    let mut output = open_for_write(outpath)?;
    for file in inputs.iter() {
        let family_name = file
            .file_stem()
            .expect("file_stem")
            .to_string_lossy()
            .to_string();

        if !families.contains(&family_name) {
            continue;
        }

        let mut reader = parse_reader(open(file)?)?;
        let basename = Path::new(file)
            .file_stem()
            .unwrap()
            .to_string_lossy()
            .to_string();

        while let Some(rec) = reader.iter_record()? {
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
    instances_dir: &PathBuf,
) -> Result<HashMap<String, PathBuf>> {
    debug!("Checking consensi_file {}", consensi_file.display());
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

    let mut consensi_names: Vec<String> = consensi_names.keys().cloned().collect();
    consensi_names.sort();
    debug!("consensi_names {consensi_names:#?}");

    let mut family_to_path: HashMap<String, PathBuf> = HashMap::new();
    for entry in fs::read_dir(instances_dir)? {
        let entry = entry?;
        let path = entry.path();
        match path.file_stem() {
            Some(stem) => {
                family_to_path.insert(stem.to_string_lossy().to_string(), path.clone())
            }
            _ => bail!("Cannot get filename from {}", path.display()),
        };
    }
    debug!("family_to_path {family_to_path:#?}");

    let instance_names: HashSet<String> = family_to_path.keys().cloned().collect();
    debug!("instance_names {instance_names:#?}");

    let missing_instances: Vec<_> = consensi_names
        .into_iter()
        .filter(|name| !instance_names.contains(name))
        .collect();

    if !missing_instances.is_empty() {
        bail!(
            "Missing the following instances: {}",
            missing_instances.iter().join(", ")
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
        consensi_names
            .insert(rec.head().trim().to_string(), rec.des().trim().to_string());
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

    // BLAST output fails to include headers
    let mut align_reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(false)
        .from_path(alignment_file)?;

    let records = align_reader.records();
    for res in records {
        let res = res?;
        let rec: RmBlastOutput = res.deserialize(None)?;
        match consensi_names.get(&rec.query) {
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
        ret.push(format!("{hours} hour{}", if hours == 1 { "" } else { "s" }));
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
    let num_from1 =
        (num_seqs_total as f64 * (num_fams1 / num_fams_total)).round() as usize;
    let num_from2 =
        (num_seqs_total as f64 * (num_fams2 / num_fams_total)).round() as usize;
    let f1 = fams1.join("::");
    let f2 = fams2.join("::");
    let new_family_name = format!("{f1}::{f2}");
    let new_family_path = args.instances_dir.join(format!("{new_family_name}.fa"));

    debug!(
        "Merging {num_from1} from {f1}, {num_from2} from {f2} => {}",
        new_family_path.display()
    );

    // Block to isolate "output" and force close when passing out of scope
    {
        let mut output = open_for_write(&new_family_path)?;
        let mut total_taken = 0;
        for (fam, num) in &[(f1, num_from1), (f2, num_from2)] {
            let fasta = family_to_path
                .get(fam)
                .unwrap_or_else(|| panic!("Missing instances for family '{fam}'"));
            total_taken += downsample(fasta, *num, num_seqs_total, &mut output)?;
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
        _ => which("Refiner").map_err(|e| anyhow!("Refiner: {e}"))?,
    };

    let mut refiner_args = vec![
        "-debug".to_string(),
        "-threads".to_string(),
        args.num_threads.to_string(),
    ];

    let rmblast = which("rmblastn").map_err(|e| anyhow!("rmblastn: {e}"))?;
    if let Some(rmblast_dir) = rmblast
        .as_path()
        .parent()
        .map(|path| path.to_string_lossy().to_string())
    {
        refiner_args.extend_from_slice(&["--rmblast_dir".to_string(), rmblast_dir]);
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

    debug!(
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
#[cfg(test)]
mod tests {
    use super::{
        bitscore_to_confidence, call_winners, cat_sequences, downsample,
        extract_scores, find_independence, format_seconds, independence,
        number_consensi, open, parse_newick, Independence, StringPair, Winners,
    };
    use anyhow::Result;
    use kseq::parse_reader;
    use pretty_assertions::assert_eq;
    use std::{
        collections::HashMap,
        fs::{self, File},
        path::PathBuf,
    };
    use tempfile::{tempdir, NamedTempFile};

    #[test]
    fn test_bitscore_to_confidence() -> Result<()> {
        let res = bitscore_to_confidence(&[&2274, &2234, &2224, &2245, &2295], 0.1227);
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
        let expected_winning_sets: HashMap<StringPair, u32> = HashMap::from([
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
        let res = cat_sequences(
            &inputs,
            &[
                "AluY".to_string(),
                "AluYa5".to_string(),
                "AluYb8".to_string(),
                "AluYb9".to_string(),
                "AluYm1".to_string(),
            ],
            &outpath,
        );
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

    //#[test]
    //fn test_check_family_instances() -> Result<()> {
    //    // The consensi contains duplicated IDs
    //    let consensi = PathBuf::from("tests/inputs/dup_consensi.fa");
    //    let instances_100 = &[PathBuf::from(
    //        "tests/inputs/instances_100/AluY.fa".to_string(),
    //    )];
    //    let res = check_family_instances(&consensi, instances_100);
    //    assert!(res.is_err());
    //    assert_eq!(
    //        res.unwrap_err().to_string(),
    //        "The following consensi IDs are duplicated: AluYb9, AluYm1"
    //    );
    //
    //    // The consensi has 5 unique IDs but the instances only 2
    //    let consensi = PathBuf::from("tests/inputs/consensi.fa");
    //    let instances_100 = &[
    //        PathBuf::from("tests/inputs/AluY.fa".to_string()),
    //        PathBuf::from("tests/inputs/AluYm1.fa".to_string()),
    //    ];
    //    let res = check_family_instances(&consensi, instances_100);
    //    assert!(res.is_err());
    //    assert_eq!(
    //        res.unwrap_err().to_string(),
    //        "Missing the following instances: AluYa5, AluYb8, AluYb9",
    //    );
    //
    //    // The consensi has 2 unique IDs but the instances have 4
    //    let consensi = PathBuf::from("tests/inputs/two_consensi.fa");
    //    let instances_100 = vec![
    //        PathBuf::from("tests/inputs/AluY.fa".to_string()),
    //        PathBuf::from("tests/inputs/AluYa5.fa".to_string()),
    //        PathBuf::from("tests/inputs/AluYb9.fa".to_string()),
    //        PathBuf::from("tests/inputs/AluYm1.fa".to_string()),
    //    ];
    //    let res = check_family_instances(&consensi, &instances_100);
    //    assert!(res.is_err());
    //    assert_eq!(
    //        res.unwrap_err().to_string(),
    //        "Missing the following consensi: AluYb9, AluYm1",
    //    );
    //
    //    Ok(())
    //}

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
        let consensi = PathBuf::from("tests/inputs/numbered_consensi.fa".to_string());
        let alignment_file = PathBuf::from("tests/inputs/blast.out");
        let prev_scores: Option<PathBuf> = None;
        let res =
            extract_scores(&alignment_file, &prev_scores, &consensi, outdir.path());
        assert!(res.is_ok());

        let scores_file = res.unwrap();
        assert!(scores_file.exists());

        let actual = fs::read_to_string(scores_file)?;
        let expected = fs::read_to_string("tests/inputs/alignment-scores.tsv")?;
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
        assert_eq!(format_seconds((one_hour * 4) + 59), "4 hours, 59 seconds");
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
        let expected = fs::read_to_string("tests/outputs/numbered_consensi.fa")?;
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
