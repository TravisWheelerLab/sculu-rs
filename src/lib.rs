use anyhow::{anyhow, bail, Result};
use assert_cmd::Command;
use clap::Parser;
use csv::{ReaderBuilder, WriterBuilder};
use itertools::Itertools;
use kseq::parse_reader;
use log::debug;
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
//use tempfile::NamedTempFile;

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

    /// PERL5LIB, e.g., to find RepeatMasker/RepeatModeler
    #[arg(long, value_name = "PERL5LIB")]
    pub perl5lib: Option<String>,

    /// Lambda value
    #[arg(long, value_name = "LAMBDA", default_value = "0.1227")]
    pub lambda: f64,

    /// Path to RepeatModeler/util/align.pl
    #[arg(long, value_name = "ALIGNER")]
    pub aligner: String,

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

    /// Debug mode
    #[arg(short, long)]
    pub debug: bool,
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

// TODO: Use type alias?
//type ClearWinnerSet = HashMap<String, u32>;
//type UnclearWinnerSet = HashMap<String, HashMap<String, u32>>;

// --------------------------------------------------
pub fn run(args: Args) -> Result<()> {
    let start = Instant::now();
    env_logger::Builder::new()
        .filter_level(if args.debug {
            log::LevelFilter::Debug
        } else {
            log::LevelFilter::Off
        })
        .init();

    debug!("args = {args:?}");
    if !args.outdir.is_dir() {
        fs::create_dir_all(&args.outdir)?;
    }

    // TODO: Check that the "consensus" sequences have clear
    // relationships to "instances" files to create "families"

    // Take the longest 100 of the instances
    let instances_100 = take_longest_100(&args)?;
    dbg!(&instances_100);

    // Align these to the consensus sequences
    let alignment_file = align(&args, &instances_100)?;

    // Extract the scores from the alignment file
    //let scores_file = PathBuf::from("sculu-out/alignment-scores.tsv");
    let scores_file = extract_scores(&alignment_file, &args)?;

    // Find the winners
    let (clear_winners, winning_sets) = call_winners(&scores_file, &args)?;

    // Find the independence of the families
    // TODO: Only return where val < 1?
    let ind_vals = independence(clear_winners, winning_sets)?;

    // Independent pairs will have a "val" of 1
    // Find the non-independent pairs
    let nonindependent: Vec<_> =
        ind_vals.iter().filter(|v| v.val < 1.).collect();

    // Merge the least independent families
    if let Some(least_ind) = nonindependent.first() {
        println!(
            "Merge {}/{} ({})",
            least_ind.f1, least_ind.f2, least_ind.val
        );
        let res = merge_families(&least_ind.f1, &least_ind.f2, &args)?;
        dbg!(res);
    }

    println!(
        r#"Finished in {} seconds, see output in "{}""#,
        start.elapsed().as_secs(),
        args.outdir.display()
    );

    Ok(())
}

// --------------------------------------------------
fn align(args: &Args, instances_100: &Vec<PathBuf>) -> Result<PathBuf> {
    // Concatente all input sequences to one file
    let all_seqs_path = args.outdir.join("all_seqs.fa");
    let mut all_seqs = open_for_write(&all_seqs_path)?;

    for filename in instances_100.iter() {
        let mut reader = parse_reader(open(filename)?)?;
        let basename = Path::new(filename)
            .file_stem()
            .unwrap()
            .to_string_lossy()
            .to_string();

        while let Some(rec) = reader.iter_record()? {
            // TODO: Remove prefix of basename?
            writeln!(
                all_seqs,
                ">{}__{}{}\n{}",
                basename,
                rec.head(),
                rec.des(),
                rec.seq()
            )?;
        }
    }

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

    debug!("Running aligner {} with {align_args:?}", args.aligner);

    let mut cmd = Command::new(&args.aligner);
    if let Some(perl5lib) = &args.perl5lib {
        cmd.env("PERL5LIB", perl5lib);
    }

    let res = cmd.args(align_args).output()?;
    if !res.status.success() {
        bail!(String::from_utf8(res.stderr)?);
    }

    let alignment_outfile = args.outdir.join("alignment.txt");
    //let mut output = File::create(&alignment_outfile)?;
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
fn extract_scores(alignment_file: &PathBuf, args: &Args) -> Result<PathBuf> {
    debug!(
        r#"Extracting scores from alignment "{}""#,
        alignment_file.display()
    );

    let file = BufReader::new(
        File::open(&alignment_file)
            .map_err(|e| anyhow!("{}: {e}", &alignment_file.display()))?,
    );

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
//) -> Result<(ClearWinnerSet, UnclearWinnerSet)> {
fn call_winners(
    scores_file: &PathBuf,
    args: &Args,
) -> Result<(HashMap<String, u32>, HashMap<String, HashMap<String, u32>>)> {
    debug!(
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
    let mut winning_sets: HashMap<String, HashMap<String, u32>> =
        HashMap::new();
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
                        let w1 = winning_sets
                            .entry(f1.to_string())
                            .or_insert(HashMap::new());
                        w1.entry(f2.to_string())
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

    Ok((clear_winners, winning_sets))
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
fn independence(
    clear_winners: HashMap<String, u32>,
    winning_sets: HashMap<String, HashMap<String, u32>>,
) -> Result<Vec<Independence>> {
    let mut families: HashSet<&str> = HashSet::new();
    for (f1, h) in winning_sets.iter() {
        families.insert(f1);
        for f2 in h.keys() {
            families.insert(f2);
        }
    }

    let mut vals = vec![];
    for pair in families.into_iter().permutations(2) {
        match pair[..] {
            [f1, f2] => {
                let num_shared = winning_sets
                    .get(f1)
                    .map_or(0u32, |v| *v.get(f2).unwrap_or(&0u32));

                let ind = if num_shared > 0 {
                    let num_wins =
                        *clear_winners.get(f1).unwrap_or(&0) as f64;
                    num_wins / (num_wins + num_shared as f64)
                } else {
                    1.
                };

                vals.push(Independence {
                    f1: f1.to_string(),
                    f2: f2.to_string(),
                    val: ind,
                });
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
fn merge_families(family1: &str, family2: &str, _args: &Args) -> Result<()> {
    debug!("Merging {family1}/{family2}");
    //let file1 = open(args.outdir.join(
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
fn take_longest_100(args: &Args) -> Result<Vec<PathBuf>> {
    debug!("Taking longest 100");

    // TODO: Make a parameter?
    let limit = 100;

    // Put subsampled sequences into outdir
    let outdir = args.outdir.join("instances_100");
    fs::create_dir_all(&outdir)?;

    let mut ret = vec![];
    for file in &args.instances {
        let min_length = find_min_len(&file, 100)?;
        let outfile = outdir.join(&file.file_name().unwrap());
        let mut output = open_for_write(&outfile)?;
        //let mut output = BufWriter::new(
        //    File::create(&outfile)
        //        .map_err(|e| anyhow!("{}: {e}", outfile.display()))?,
        //);
        let mut reader = parse_reader(open(&file)?)?;
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
#[cfg(test)]
mod tests {
    use super::{
        align, bitscore_to_confidence, call_winners, extract_scores,
        independence, Args,
    };
    use anyhow::Result;
    use pretty_assertions::assert_eq;
    use std::{collections::HashMap, fs, path::PathBuf};
    use tempfile::tempdir;

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
            align_gap_init: -25,
            align_extension: -5,
            align_min_match: 7,
            align_bandwidth: 14,
            align_mask_level: 101,
            align_min_score: 200,
            debug: false,
        };
        let instances = vec![
            PathBuf::from("tests/inputs/instances_100/AluY.fa"),
            PathBuf::from("tests/inputs/instances_100/AluYa5.fa"),
            PathBuf::from("tests/inputs/instances_100/AluYb8.fa"),
            PathBuf::from("tests/inputs/instances_100/AluYb9.fa"),
            PathBuf::from("tests/inputs/instances_100/AluYm1.fa"),
        ];
        let res = align(&args, &instances);
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
            align_gap_init: -25,
            align_extension: -5,
            align_min_match: 7,
            align_bandwidth: 14,
            align_mask_level: 101,
            align_min_score: 200,
            debug: false,
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
            align_gap_init: -25,
            align_extension: -5,
            align_min_match: 7,
            align_bandwidth: 14,
            align_mask_level: 101,
            align_min_score: 200,
            debug: false,
        };
        let res = call_winners(&scores_file, &args);
        assert!(res.is_ok());

        let (clear_winners, winning_sets) = res.unwrap();

        assert_eq!(
            clear_winners,
            HashMap::from([
                ("AluY".to_string(), 37),
                ("AluYm1".to_string(), 36),
                ("AluYb8".to_string(), 50)
            ])
        );

        assert_eq!(
            winning_sets,
            HashMap::from([
                (
                    "AluYa5".to_string(),
                    HashMap::from([
                        ("AluYb8".to_string(), 6u32),
                        ("AluYm1".to_string(), 4u32),
                    ]),
                ),
                (
                    "AluY".to_string(),
                    HashMap::from([
                        ("AluYb8".to_string(), 6),
                        ("AluYa5".to_string(), 12),
                        ("AluYm1".to_string(), 46),
                    ]),
                ),
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
            (
                "AluYa5".to_string(),
                HashMap::from([
                    ("AluYb8".to_string(), 6u32),
                    ("AluYm1".to_string(), 4u32),
                ]),
            ),
            (
                "AluY".to_string(),
                HashMap::from([
                    ("AluYb8".to_string(), 6),
                    ("AluYa5".to_string(), 12),
                    ("AluYm1".to_string(), 46),
                ]),
            ),
        ]);

        let res = independence(clear_winners, winning_sets);
        assert!(res.is_ok());
        Ok(())
    }
}
