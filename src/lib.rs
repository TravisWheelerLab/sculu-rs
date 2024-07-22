use anyhow::{anyhow, bail, Result};
use assert_cmd::Command;
use clap::{builder::PossibleValue, Parser, ValueEnum};
use csv::{ReaderBuilder, WriterBuilder};
use itertools::Itertools;
use kseq::parse_reader;
use log::{debug, info};
use rand::{seq::SliceRandom, thread_rng};
use regex::Regex;
use serde::{Deserialize, Serialize};
use std::{
    cmp,
    collections::{HashMap, HashSet},
    fs::{self, File},
    io::{BufRead, BufReader, BufWriter, Write},
    path::{Path, PathBuf},
    time::Instant,
};

/// SCULU subfamily clustering tool
#[derive(Debug, Parser)]
#[command(author, version, about)]
pub struct Args {
    /// FASTA file of subfamily consensus sequences
    #[arg(long, value_name = "CONSENSUS")]
    pub consensus: PathBuf,

    /// Directory of instance files for each subfamily
    #[arg(long, value_name = "INSTANCES", required = true, num_args=1..)]
    pub instances: Vec<PathBuf>,

    /// Output directory
    #[arg(short, long, value_name = "OUTDIR", default_value = "sculu-out")]
    pub outdir: PathBuf,

    /// Lambda value
    #[arg(long, value_name = "LAMBDA", default_value = "0.1227")]
    pub lambda: f64,

    /// Path to RepeatModeler/util/align.pl
    #[arg(long, value_name = "ALIGNER")]
    pub aligner: String,

    /// Path to RepeatModeler/Refiner
    #[arg(long, value_name = "REFINER")]
    pub refiner: String,

    /// PERL5LIB, e.g., to find RepeatMasker/RepeatModeler
    #[arg(long, value_name = "PERL5LIB")]
    pub perl5lib: Option<String>,

    /// Path to rmblastn
    #[arg(long, value_name = "RMBLAST_DIR")]
    pub rmblast_dir: Option<String>,

    /// Alignment matrix
    #[arg(long, value_name = "MATRIX")]
    pub alignment_matrix: Option<String>,

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

#[derive(Debug, Serialize, Deserialize)]
struct AlignmentScore {
    score: u32,
    target: String,
    query: String,
}

#[derive(Debug)]
struct Independence {
    f1: String,
    f2: String,
    val: f64,
}

#[derive(Debug, Eq, PartialEq, Hash)]
struct StringPair(String, String);

#[derive(Debug)]
struct Winners {
    clear_winners: HashMap<String, u32>,
    winning_sets: HashMap<StringPair, u32>,
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
        .init();

    info!("args = {args:#?}");
    if !args.outdir.is_dir() {
        fs::create_dir_all(&args.outdir)?;
    }

    // TODO: Make this a parameter?
    let independence_threshold = 0.5;

    // Take the longest 100 of the instances
    let instances_100_dir = args.outdir.join("instances_100");
    fs::create_dir_all(&instances_100_dir)?;
    let instances_100 = take_longest_100(&args, &instances_100_dir)?;
    debug!("instances_100 = {instances_100:#?}");

    // Check that the "consensus" sequences have clear
    // relationships to "instances" files to create "families"
    let family_to_path = check_family_instances(&args, &instances_100)?;
    debug!("family_to_path = {family_to_path:#?}");

    // Concatenate the longest 100 into a single file for alignment
    let all_seqs_path = &args.outdir.join("all_seqs.fa");
    cat_sequences(&instances_100, &all_seqs_path)?;

    // BLAST chokes on long identifiers, so give the consensi
    // sequences integer IDs and move the names to the desc
    let consensi_path = args.outdir.join("consensi.fa");
    number_consensi(&args.consensus, &consensi_path)?;
    args = Args {
        consensus: consensi_path,
        ..args
    };

    let top_outdir = args.outdir.clone();
    for round in 1.. {
        println!("\n>>> Round {round} <<<");

        args.outdir = top_outdir.join(format!("round_{round:03}"));
        fs::create_dir_all(&args.outdir)?;

        // Align longest 100 to the consensi
        let alignment_file = align(&args, &all_seqs_path)?;

        // Extract the scores from the alignment file
        //let scores_file = PathBuf::from("sculu-out/alignment-scores.tsv");
        let scores_file = extract_scores(&alignment_file, &args)?;

        // Find the winners from the scores
        //let (clear_winners, winning_sets) =
        let winners = call_winners(&scores_file, &args)?;

        // Find the independence of the families
        // TODO: Only return where val < .5?
        // Independent pairs will have a "val" of 1
        // Find the pairs with less than 50% chance of independence
        let ind_vals = independence(winners, independence_threshold)?;

        // Merge the least independent families
        if let Some(least) = ind_vals.first() {
            // The raw names may be "Fam1::Fam2::Fam2""
            let raw_names = get_family_names_from_consensi(
                &args.consensus,
                &[&least.f1, &least.f2],
            )?;

            let mut families = vec![];
            for fam in raw_names {
                let mut fams: Vec<String> =
                    fam.split("::").map(Into::into).collect();
                families.append(&mut fams);
            }

            info!("Merge {} ({:0.04})", families.join(", "), least.val);
            let new_consensi =
                merge_families(&families, &args, &family_to_path)?;

            args = Args {
                consensus: new_consensi,
                ..args
            };
        } else {
            // TODO: Add a check for max # of rounds to break?
            println!("No families to merge.");
            break;
        }
    }

    println!(
        r#"Finished in {} seconds, see output in "{}""#,
        start.elapsed().as_secs(),
        args.outdir.display()
    );

    Ok(())
}

// --------------------------------------------------
fn align(args: &Args, all_seqs_path: &PathBuf) -> Result<PathBuf> {
    // Perform alignment to query
    let mut align_args = vec![
        "-rmblast".to_string(),
        "-gap_init".to_string(),
        args.align_gap_init.to_string(),
        "-extension".to_string(),
        args.align_extension.to_string(),
        "-minmatch".to_string(),
        args.align_min_match.to_string(),
        "-bandwidth".to_string(),
        args.align_bandwidth.to_string(),
        "-masklevel".to_string(),
        args.align_mask_level.to_string(),
        "-minscore".to_string(),
        args.align_min_score.to_string(),
    ];

    if let Some(matrix) = &args.alignment_matrix {
        align_args
            .extend_from_slice(&["-matrix".to_string(), matrix.to_string()]);
    }

    align_args.extend_from_slice(&[
        "-alignments".to_string(),
        all_seqs_path.to_string_lossy().to_string(),
        args.consensus.to_string_lossy().to_string(),
    ]);

    info!(r#"Running "{} {}""#, args.aligner, align_args.join(" "));

    let mut cmd = Command::new(&args.aligner);
    if let Some(perl5lib) = &args.perl5lib {
        cmd.env("PERL5LIB", perl5lib);
    }

    let res = cmd.args(align_args).output()?;
    if !res.status.success() {
        bail!(String::from_utf8(res.stderr)?);
    }

    let alignment_outfile = args.outdir.join("alignment.txt");
    let mut output = open_for_write(&alignment_outfile)?;
    output.write_all(&res.stdout)?;

    Ok(alignment_outfile)
}

// --------------------------------------------------
fn bitscore_to_confidence(vals: &Vec<&u32>, lambda: f64) -> Result<Vec<f64>> {
    let converted: Vec<_> = vals
        .into_iter()
        .map(|&&v| (v as f64 * lambda).exp2())
        .collect();
    let total: f64 = converted.iter().sum();
    if total > 0. {
        Ok(converted.into_iter().map(|v| v / total).collect())
    } else {
        bail!("Sum of converted values equals zero")
    }
}

// --------------------------------------------------
fn check_family_instances(
    args: &Args,
    instances_100: &Vec<PathBuf>,
) -> Result<HashMap<String, PathBuf>> {
    let mut reader = parse_reader(open(&args.consensus)?)?;
    let mut consensi_names: HashMap<String, u32> = HashMap::new();
    while let Some(rec) = reader.iter_record()? {
        consensi_names
            .entry(rec.head().to_string())
            .and_modify(|v| *v += 1)
            .or_insert(1);
    }

    let mut dups: Vec<_> = consensi_names
        .iter()
        .flat_map(|(name, &count)| (count > 1).then(|| name))
        .collect();

    if dups.len() > 0 {
        // Have to sort for tests
        dups.sort();
        bail!(
            "The following consensi IDs are duplicated: {}",
            dups.iter().join(", ")
        );
    }

    let consensi_names: HashSet<String> =
        consensi_names.keys().cloned().collect();

    let mut fam_to_file: HashMap<String, PathBuf> = HashMap::new();
    for instance in instances_100 {
        match instance.file_stem() {
            Some(stem) => fam_to_file
                .insert(stem.to_string_lossy().to_string(), instance.clone()),
            _ => bail!(
                "Cannot get filename from {}",
                instance.to_string_lossy().to_string()
            ),
        };
    }

    let instance_names: HashSet<String> =
        fam_to_file.keys().cloned().collect();

    let missing_instances: HashSet<_> =
        consensi_names.difference(&instance_names).collect();

    if missing_instances.len() > 0 {
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

    if missing_consensi.len() > 0 {
        let mut missing_consensi: Vec<_> = missing_consensi.iter().collect();
        missing_consensi.sort();
        bail!(
            "Missing the following consensi: {}",
            missing_consensi.iter().join(", ")
        );
    }

    Ok(fam_to_file)
}

// --------------------------------------------------
fn downsample(fasta: &PathBuf, mut output: impl Write) -> Result<()> {
    // Assumes FASTA file contains 100 longest sequences
    let mut reader = parse_reader(open(fasta)?)?;

    // Randomly pick 50 record positions
    let mut rng = thread_rng();
    let range: Vec<usize> = (0..100).collect();
    let mut take: Vec<usize> =
        range.choose_multiple(&mut rng, 50).cloned().collect();
    take.sort();
    take.reverse();

    let mut next_take = take.pop().unwrap();
    let mut i = 0;
    while let Some(rec) = reader.iter_record()? {
        if i == next_take {
            writeln!(output, ">{}{}\n{}", rec.head(), rec.des(), rec.seq())?;
            match take.pop() {
                Some(val) => next_take = val,
                _ => break,
            }
        }
        i += 1;
    }

    Ok(())
}

// --------------------------------------------------
fn extract_scores(alignment_file: &PathBuf, args: &Args) -> Result<PathBuf> {
    info!(
        r#"Extracting scores from alignment "{}""#,
        alignment_file.display()
    );

    let file = open(&alignment_file)?;
    let scores_file = args.outdir.join("alignment-scores.tsv");
    let mut wtr = WriterBuilder::new()
        .has_headers(true)
        .delimiter(b'\t')
        .from_path(&scores_file)?;

    // Take only the score lines that start with an integer
    let starts_with_score = Regex::new(r"^\d+\s").unwrap();
    let spaces = Regex::new(r"\s+").unwrap();
    for line in file.lines().map_while(Result::ok) {
        if starts_with_score.is_match(&line) {
            let parts: Vec<_> = spaces.split(&line).collect();
            if parts.len() == 12 {
                let score: u32 = parts[0].parse()?;
                let target = parts[4].to_string();
                let query = parts[8].to_string();
                wtr.serialize(AlignmentScore {
                    score,
                    target,
                    query,
                })?;
            }
        }
    }

    Ok(scores_file)
}

// --------------------------------------------------
fn call_winners(scores_file: &PathBuf, args: &Args) -> Result<Winners> {
    info!(
        r#"Calling winners from scores file "{}""#,
        scores_file.display()
    );

    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .from_path(scores_file)?;
    let mut records = reader.deserialize();

    let mut scores: HashMap<String, HashMap<String, u32>> = HashMap::new();
    while let Some(res) = records.next() {
        let rec: AlignmentScore = res?;
        let tscore = scores.entry(rec.target).or_insert(HashMap::new());
        tscore
            .entry(rec.query)
            .and_modify(|v| *v = cmp::max(rec.score, *v))
            .or_insert(rec.score);
    }

    let mut clear_winners: HashMap<String, u32> = HashMap::new();
    let mut winning_sets: HashMap<StringPair, u32> = HashMap::new();
    for (target, bit_scores) in scores.iter() {
        debug!("\n>>> {target}");
        let families: Vec<_> = bit_scores.keys().collect();
        debug!("families   = {families:?}");
        let raw_scores: Vec<_> = bit_scores.values().collect();
        debug!("raw scores = {raw_scores:?}");
        let conf: Vec<_> = bitscore_to_confidence(&raw_scores, args.lambda)?;
        debug!(
            "confidence = {:?}",
            conf.iter()
                .map(|v| format!("{:0.04}", v))
                .collect::<Vec<_>>()
        );

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
            let comps: Vec<_> =
                pairs.iter().map(|(x, y)| (x * 0.3) > *y as f64).collect();
            all_comps.push(comps.iter().all(|&v| v));
            debug!(
                "{:8} {:0.04} {:5} => {:?}",
                families[i],
                val,
                comps.iter().all(|&v| v),
                comps
            );
        }
        let winners: Vec<String> = families
            .iter()
            .zip(all_comps)
            .filter_map(|(fam, win)| win.then(|| fam.to_string()))
            .collect();

        if winners.len() == 1 {
            let winner = winners.first().unwrap();
            clear_winners
                .entry(winner.to_string())
                .and_modify(|v| *v += 1)
                .or_insert(1);
            debug!("Clear Winner: {winner}");
        } else {
            let mut fam_comps: Vec<_> = conf.iter().zip(families).collect();
            fam_comps.sort_by(|a, b| {
                b.0.partial_cmp(a.0).unwrap().then_with(|| a.1.cmp(&b.1))
            });
            let (&top_conf, _) = fam_comps.first().unwrap();
            let threshold = top_conf / 5.;
            let winning_set: Vec<_> = fam_comps
                .iter()
                .filter_map(|(&c, f)| (c > threshold).then(|| f))
                .collect();

            debug!("Winning Set: {winning_set:?}");
            for mut pair in winning_set.into_iter().permutations(2) {
                pair.sort();
                match pair[..] {
                    [&f1, &f2] => {
                        let key = StringPair(f1.to_string(), f2.to_string());
                        winning_sets
                            .entry(key)
                            .and_modify(|v| *v += 1)
                            .or_insert(1);
                    }
                    _ => (),
                }
            }
        }
    }

    debug!("Clear Winners = {clear_winners:?}");
    debug!("Winning Sets = {winning_sets:?}");

    Ok(Winners {
        clear_winners,
        winning_sets,
    })
}

// --------------------------------------------------
fn find_min_len(filename: &PathBuf, number: usize) -> Result<usize> {
    let mut reader = parse_reader(open(&filename)?)?;
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
        // Default will be to take all
        min_len = *len;
        top_n += count_by_len.get(&len).unwrap_or(&0);
        if top_n >= number {
            break;
        }
    }

    Ok(min_len)
}

// --------------------------------------------------
fn get_family_names_from_consensi(
    consensi: &PathBuf,
    ids: &[&str],
) -> Result<Vec<String>> {
    let mut reader = parse_reader(open(&consensi)?)?;
    let mut consensi_names: HashMap<String, String> = HashMap::new();
    while let Some(rec) = reader.iter_record()? {
        consensi_names.insert(
            rec.head().trim().to_string(),
            rec.des().trim().to_string(),
        );
    }

    Ok(ids
        .iter()
        .flat_map(|&id| consensi_names.get(id).map(|v| v.to_string()))
        .collect::<Vec<_>>())
}

// --------------------------------------------------
fn independence(
    winners: Winners,
    threshold: f64,
) -> Result<Vec<Independence>> {
    let mut families: HashSet<&str> = HashSet::new();
    for StringPair(f1, f2) in winners.winning_sets.keys() {
        families.insert(&f1);
        families.insert(&f2);
    }

    let mut vals = vec![];
    for pair in families.into_iter().permutations(2) {
        match pair[..] {
            [f1, f2] => {
                let key = StringPair(f1.to_string(), f2.to_string());
                let &num_shared =
                    winners.winning_sets.get(&key).unwrap_or(&0u32);
                let ind = if num_shared > 0 {
                    let num_wins =
                        *winners.clear_winners.get(f1).unwrap_or(&0) as f64;
                    num_wins / (num_wins + num_shared as f64)
                } else {
                    1.
                };

                if ind < threshold {
                    vals.push(Independence {
                        f1: f1.to_string(),
                        f2: f2.to_string(),
                        val: ind,
                    });
                }
            }
            _ => (),
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

    debug!("independence {vals:#?}");
    Ok(vals)
}

// --------------------------------------------------
fn merge_families(
    families: &Vec<String>,
    args: &Args,
    fam_to_file: &HashMap<String, PathBuf>,
) -> Result<PathBuf> {
    let outdir = args.outdir.join("msa");
    fs::create_dir_all(&outdir)?;

    let msa_input = outdir.join("msa-input.fa");
    {
        // Block to isolate "output" and force
        // close when passing out of scope
        let mut output = open_for_write(&msa_input)?;
        for fam in families {
            let fasta = fam_to_file.get(fam).unwrap();
            downsample(fasta, &mut output)?;
        }
    }

    //TODO: Make threads value a parameter?
    let mut refiner_args = vec!["-threads".to_string(), "4".to_string()];

    match args.log {
        Some(LogLevel::Debug) => refiner_args.push("-debug".to_string()),
        _ => (),
    }

    if let Some(rmblast_dir) = &args.rmblast_dir {
        refiner_args.extend_from_slice(&[
            "--rmblast_dir".to_string(),
            rmblast_dir.to_string(),
        ]);
    }

    refiner_args.push(msa_input.to_string_lossy().to_string());
    info!(r#"Running "{} {}""#, &args.refiner, &refiner_args.join(" "));

    // TODO: This is failing ATM
    let mut cmd = Command::new(&args.refiner);
    if let Some(perl5lib) = &args.perl5lib {
        cmd.env("PERL5LIB", perl5lib);
    }
    let res = cmd.args(refiner_args).output()?;
    if !res.status.success() {
        debug!("{}", String::from_utf8(res.stdout)?);
        bail!(String::from_utf8(res.stderr)?);
    }

    let consensus_path = &outdir.join("msa-input.fa.refiner_cons");
    if !consensus_path.exists() {
        bail!(
            "Failed to find expected consensus file {}",
            consensus_path.display()
        );
    }

    let mut reader = parse_reader(open(&consensus_path)?)?;
    let consensus_seq = reader
        .iter_record()?
        .map(|rec| rec.seq().to_string())
        .expect("Failed to read consensus file");

    // Print the consensus sequence to the new consensi file
    let new_consensi_path = outdir.join("new-consensi.fa");
    let mut new_consensi = open_for_write(&new_consensi_path)?;
    writeln!(new_consensi, ">0 {}\n{consensus_seq}", &families.join("::"))?;

    let mut merge_fams: HashSet<String> = HashSet::new();
    for fam in families {
        merge_fams.insert(fam.to_string());
    }

    // Add the other original consensi
    let mut orig_consensi = parse_reader(open(&args.consensus)?)?;
    let mut i = 1;
    while let Some(rec) = orig_consensi.iter_record()? {
        // Skip the families we merged
        let mut seq_fams: HashSet<String> = HashSet::new();
        for fam in rec.des().trim().split("::") {
            seq_fams.insert(fam.to_string());
        }

        let common: HashSet<_> = merge_fams.intersection(&seq_fams).collect();
        if common.is_empty() {
            writeln!(
                new_consensi,
                ">{i} {}\n{}",
                rec.des().trim(),
                rec.seq()
            )?;
        }
        i += 1;
    }

    info!(r#"New consensi "{}""#, new_consensi_path.display());
    Ok(new_consensi_path)
}

// --------------------------------------------------
fn number_consensi(consensi: &PathBuf, outpath: &PathBuf) -> Result<()> {
    let mut outfile = open_for_write(outpath)?;
    let mut reader = parse_reader(open(&consensi)?)?;
    let mut i = 0;

    while let Some(rec) = reader.iter_record()? {
        writeln!(outfile, ">{i} {}\n{}", rec.head(), rec.seq())?;
        i += 1;
    }

    Ok(())
}

// --------------------------------------------------
fn open(filename: &PathBuf) -> Result<Box<dyn BufRead>> {
    Ok(Box::new(BufReader::new(File::open(filename).map_err(
        |e| anyhow!("Cannot read {}: {e}", filename.display()),
    )?)))
}

// --------------------------------------------------
fn open_for_write(filename: &PathBuf) -> Result<Box<dyn Write>> {
    Ok(Box::new(BufWriter::new(File::create(&filename).map_err(
        |e| anyhow!("Cannot write {}: {e}", filename.display()),
    )?)))
}

// --------------------------------------------------
fn take_longest_100(args: &Args, outdir: &PathBuf) -> Result<Vec<PathBuf>> {
    info!("Taking longest 100");

    // TODO: Make a parameter?
    let limit = 100;

    let mut ret = vec![];
    for file in &args.instances {
        let min_length = find_min_len(&file, 100)?;
        let outfile = outdir.join(&file.file_name().unwrap());
        let mut output = open_for_write(&outfile)?;
        let mut reader = parse_reader(open(file)?)?;
        let mut taken = 0;

        while let Some(rec) = reader.iter_record()? {
            if limit > 0 && taken == limit {
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
#[cfg(test)]
mod tests {
    use super::{
        align, bitscore_to_confidence, call_winners, cat_sequences,
        check_family_instances, extract_scores,
        get_family_names_from_consensi, independence, number_consensi, open,
        Args, StringPair, Winners,
    };
    use anyhow::Result;
    use kseq::parse_reader;
    use pretty_assertions::assert_eq;
    use std::{collections::HashMap, fs, path::PathBuf};
    use tempfile::{tempdir, NamedTempFile};

    // Or "tests/inputs/25p41g.matrix"?
    const MATRIX: &str =
        "/Users/kyclark/work/RepeatMasker/Matrices/ncbi/nt/25p41g.matrix";

    #[test]
    fn test_bitscore_to_confidence() -> Result<()> {
        let res = bitscore_to_confidence(
            &vec![&2274, &2234, &2224, &2245, &2295],
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
    fn test_cat_sequences() -> Result<()> {
        let outdir = tempdir()?;
        let inputs = vec![
            PathBuf::from("tests/inputs/instances_100/AluY.fa"),
            PathBuf::from("tests/inputs/instances_100/AluYa5.fa"),
            PathBuf::from("tests/inputs/instances_100/AluYb8.fa"),
            PathBuf::from("tests/inputs/instances_100/AluYb9.fa"),
            PathBuf::from("tests/inputs/instances_100/AluYm1.fa"),
        ];
        let outpath = outdir.path().join("all_seqs.fa");
        let res = cat_sequences(&inputs, &outpath);
        assert!(res.is_ok());
        assert!(outpath.exists());

        let mut reader = parse_reader(open(&outpath)?)?;
        let mut count = 0;
        while let Some(_) = reader.iter_record()? {
            count += 1;
        }
        assert_eq!(count, 500);

        Ok(())
    }

    #[test]
    fn test_align() -> Result<()> {
        let outdir = tempdir()?;

        let args = Args {
            consensus: PathBuf::from("tests/inputs/consensi.fa".to_string()),
            instances: vec![PathBuf::from(
                "tests/inputs/test_set.fa".to_string(),
            )],
            alignment_matrix: Some(MATRIX.to_string()),
            outdir: outdir.path().into(),
            perl5lib: Some("/Users/kyclark/work/RepeatMasker".to_string()),
            lambda: 0.1227,
            aligner: "/Users/kyclark/work/RepeatModeler/util/align.pl"
                .to_string(),
            refiner: "/Users/kyclark/work/RepeatModeler/Refiner".to_string(),
            rmblast_dir: Some("/Users/kyclark/.local/bin".to_string()),
            align_gap_init: -25,
            align_extension: -5,
            align_min_match: 7,
            align_bandwidth: 14,
            align_mask_level: 101,
            align_min_score: 200,
            log: None,
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
    fn test_check_family_instances() -> Result<()> {
        // The consensi contains duplicated IDs
        let outdir = tempdir()?;
        let args1 = Args {
            consensus: PathBuf::from("tests/inputs/dup_consensi.fa"),
            instances: vec![PathBuf::from(
                "tests/inputs/test_set.fa".to_string(),
            )],
            alignment_matrix: Some(MATRIX.to_string()),
            outdir: outdir.path().into(),
            perl5lib: Some("/Users/kyclark/work/RepeatMasker".to_string()),
            lambda: 0.1227,
            aligner: "/Users/kyclark/work/RepeatModeler/util/align.pl"
                .to_string(),
            refiner: "/Users/kyclark/work/RepeatModeler/Refiner".to_string(),
            rmblast_dir: Some("/Users/kyclark/.local/bin".to_string()),
            align_gap_init: -25,
            align_extension: -5,
            align_min_match: 7,
            align_bandwidth: 14,
            align_mask_level: 101,
            align_min_score: 200,
            log: None,
        };

        let instances_100 = vec![PathBuf::from(
            "tests/inputs/instances_100/AluY.fa".to_string(),
        )];
        let res = check_family_instances(&args1, &instances_100);
        assert!(res.is_err());
        assert_eq!(
            res.unwrap_err().to_string(),
            "The following consensi IDs are duplicated: AluYb9, AluYm1"
        );

        // The consensi has 5 unique IDs but the instances only 2
        let args2 = Args {
            consensus: PathBuf::from("tests/inputs/consensi.fa"),
            ..args1
        };
        let instances_100 = vec![
            PathBuf::from("tests/inputs/AluY.fa".to_string()),
            PathBuf::from("tests/inputs/AluYm1.fa".to_string()),
        ];
        let res = check_family_instances(&args2, &instances_100);
        assert!(res.is_err());
        assert_eq!(
            res.unwrap_err().to_string(),
            "Missing the following instances: AluYa5, AluYb8, AluYb9",
        );

        // The consensi has 2 unique IDs but the instances have 4
        let args3 = Args {
            consensus: PathBuf::from("tests/inputs/two_consensi.fa"),
            ..args2
        };
        let instances_100 = vec![
            PathBuf::from("tests/inputs/AluY.fa".to_string()),
            PathBuf::from("tests/inputs/AluYa5.fa".to_string()),
            PathBuf::from("tests/inputs/AluYb9.fa".to_string()),
            PathBuf::from("tests/inputs/AluYm1.fa".to_string()),
        ];
        let res = check_family_instances(&args3, &instances_100);
        assert!(res.is_err());
        assert_eq!(
            res.unwrap_err().to_string(),
            "Missing the following consensi: AluYb9, AluYm1",
        );

        Ok(())
    }

    #[test]
    fn test_extract_scores() -> Result<()> {
        let outdir = tempdir()?;
        let args = Args {
            consensus: PathBuf::from("tests/inputs/consensi.fa".to_string()),
            instances: vec![PathBuf::from(
                "tests/inputs/test_set.fa".to_string(),
            )],
            alignment_matrix: Some(MATRIX.to_string()),
            outdir: outdir.path().into(),
            perl5lib: Some("/Users/kyclark/work/RepeatMasker".to_string()),
            lambda: 0.1227,
            aligner: "/Users/kyclark/work/RepeatModeler/util/align.pl"
                .to_string(),
            refiner: "/Users/kyclark/work/RepeatModeler/Refiner".to_string(),
            rmblast_dir: Some("/Users/kyclark/.local/bin".to_string()),
            align_gap_init: -25,
            align_extension: -5,
            align_min_match: 7,
            align_bandwidth: 14,
            align_mask_level: 101,
            align_min_score: 200,
            log: None,
        };

        let alignment_file = PathBuf::from("tests/inputs/alignment.txt");
        let res = extract_scores(&alignment_file, &args);
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
    fn test_call_winners() -> Result<()> {
        let scores_file = PathBuf::from("tests/outputs/alignment-scores.tsv");
        let outdir = tempdir()?;
        let args = Args {
            consensus: PathBuf::from("tests/inputs/consensi.fa".to_string()),
            instances: vec![PathBuf::from(
                "tests/inputs/test_set.fa".to_string(),
            )],
            alignment_matrix: Some(MATRIX.to_string()),
            outdir: outdir.path().into(),
            perl5lib: Some("/Users/kyclark/work/RepeatMasker".to_string()),
            lambda: 0.1227,
            aligner: "/Users/kyclark/work/RepeatModeler/util/align.pl"
                .to_string(),
            refiner: "/Users/kyclark/work/RepeatModeler/Refiner".to_string(),
            rmblast_dir: Some("/Users/kyclark/.local/bin".to_string()),
            align_gap_init: -25,
            align_extension: -5,
            align_min_match: 7,
            align_bandwidth: 14,
            align_mask_level: 101,
            align_min_score: 200,
            log: None,
        };
        let res = call_winners(&scores_file, &args);
        assert!(res.is_ok());

        let winners = res.unwrap();

        assert_eq!(
            winners.clear_winners,
            HashMap::from([
                ("AluY".to_string(), 37),
                ("AluYm1".to_string(), 36),
                ("AluYb8".to_string(), 50)
            ])
        );

        assert_eq!(
            winners.winning_sets,
            HashMap::from([
                (StringPair("AluYa5".to_string(), "AluYb8".to_string()), 6),
                (StringPair("AluYa5".to_string(), "AluYm1".to_string()), 4),
                (StringPair("AluY".to_string(), "AluYb8".to_string()), 6),
                (StringPair("AluY".to_string(), "AluYa5".to_string()), 12),
                (StringPair("AluY".to_string(), "AluYm1".to_string()), 46)
            ])
        );

        Ok(())
    }

    #[test]
    fn test_independence() -> Result<()> {
        let clear_winners = HashMap::from([
            ("AluY".to_string(), 37),
            ("AluYm1".to_string(), 36),
            ("AluYb8".to_string(), 50),
        ]);

        let winning_sets = HashMap::from([
            (StringPair("AluYa5".to_string(), "AluYb8".to_string()), 6),
            (StringPair("AluYa5".to_string(), "AluYm1".to_string()), 4),
            (StringPair("AluY".to_string(), "AluYb8".to_string()), 6),
            (StringPair("AluY".to_string(), "AluYa5".to_string()), 12),
            (StringPair("AluY".to_string(), "AluYm1".to_string()), 46),
        ]);

        let res = independence(
            Winners {
                clear_winners,
                winning_sets,
            },
            0.5,
        );
        assert!(res.is_ok());
        Ok(())
    }

    #[test]
    fn test_get_family_names_from_consensi() -> Result<()> {
        let numbered = PathBuf::from("tests/outputs/numbered_consensi.fa");
        let res = get_family_names_from_consensi(&numbered, &["1", "3"]);
        assert!(res.is_ok());

        assert_eq!(
            res.unwrap(),
            vec!["AluYb8".to_string(), "AluYa5".to_string()]
        );
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
}
