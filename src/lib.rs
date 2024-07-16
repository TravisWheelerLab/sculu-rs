use anyhow::{anyhow, bail, Result};
use assert_cmd::Command;
use clap::Parser;
use itertools::Itertools;
use kseq::parse_reader;
use regex::Regex;
use std::{
    cmp,
    collections::HashMap,
    fs::File,
    io::{BufReader, Write},
    path::Path,
};
//use tempfile::NamedTempFile;

/// SCULU subfamily clustering tool
#[derive(Debug, Parser)]
#[command(author, version, about)]
pub struct Args {
    /// FASTA file of subfamily consensus sequences
    #[arg(long, value_name = "CONSENSUS")]
    pub consensus: String,

    /// Directory of instance files for each subfamily
    #[arg(long, value_name = "INSTANCES", required = true, num_args=1..)]
    pub instances: Vec<String>,

    /// PERL5LIB, e.g., to find RepeatMasker/RepeatModeler
    #[arg(long, value_name = "PERL5LIB")]
    pub perl5lib: Option<String>,

    /// Path to RepeatModeler/util/align.pl
    #[arg(long, value_name = "ALIGNER")]
    pub aligner: String,

    /// Lambda value
    #[arg(long, value_name = "LAMBDA", default_value = "0.1227")]
    pub lambda: f64,
}

// --------------------------------------------------
pub fn run(args: Args) -> Result<()> {
    dbg!(&args);
    let _align_res = align(&args)?;
    Ok(())
}

// --------------------------------------------------
fn align(args: &Args) -> Result<()> {
    //let mut all_seqs = NamedTempFile::new()?;
    //let all_seqs_path = &all_seqs.path().to_str().unwrap();

    let all_seqs_path = "all_seqs.fa";
    let mut all_seqs = File::create(all_seqs_path)?;
    println!("Writing to {all_seqs_path}");
    for filename in &args.instances {
        let file = BufReader::new(
            File::open(filename).map_err(|e| anyhow!("{}: {e}", filename))?,
        );
        let mut reader = parse_reader(file)?;
        let basename = Path::new(filename)
            .file_stem()
            .unwrap()
            .to_string_lossy()
            .to_string();
        while let Some(rec) = reader.iter_record()? {
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

    let align_args = &[
        "-rmblast",
        "-gap_init",
        "-25",
        "-extension",
        "-5",
        "-minmatch",
        "7",
        "-bandwidth",
        "14",
        "-masklevel",
        "101",
        "-matrix",
        "/Users/kyclark/work/RepeatMasker/Matrices/ncbi/nt/25p41g.matrix",
        "-minscore",
        "200",
        "-alignments",
        all_seqs_path,
        &args.consensus,
    ];
    println!("{} {:?}", args.aligner, align_args);

    let mut cmd = Command::new(&args.aligner);
    if let Some(perl5lib) = &args.perl5lib {
        cmd.env("PERL5LIB", perl5lib);
    }

    //let mut align_out = NamedTempFile::new()?;
    //let all_seqs_path = &all_seqs.path().to_str().unwrap();
    //
    //let align_out = "align.out";
    //let out = File::create(align_out)?;
    let res = cmd.args(align_args).output()?;
    if !res.status.success() {
        bail!(String::from_utf8(res.stderr)?);
    }
    let stdout = String::from_utf8(res.stdout)?;

    let mut scores: HashMap<String, HashMap<String, u32>> = HashMap::new();
    let starts_with_score = Regex::new(r"^\d+\s").unwrap();
    let spaces = Regex::new(r"\s+").unwrap();
    for line in stdout.split('\n') {
        if starts_with_score.is_match(line) {
            let parts: Vec<_> = spaces.split(line).collect();
            if parts.len() == 12 {
                let score: u32 = parts[0].parse()?;
                let target = parts[4].to_string();
                let query = parts[8].to_string();
                // println!("target {target} query {query} score {score}");

                let tscore = scores.entry(target).or_insert(HashMap::new());
                tscore
                    .entry(query)
                    .and_modify(|v| *v = cmp::max(score, *v))
                    .or_insert(score);
            }
        }
    }
    //println!("{scores:#?}");

    let mut clear_winners: HashMap<String, u32> = HashMap::new();
    let mut winning_sets: HashMap<String, HashMap<String, u32>> =
        HashMap::new();
    for (target, bit_scores) in scores.iter() {
        println!(">>> {target}");
        let families: Vec<_> = bit_scores.keys().collect();
        println!("{families:?}");
        let conf: Vec<_> = bitscore_to_confidence(
            &bit_scores.values().collect::<Vec<_>>(),
            args.lambda,
        )?;
        println!("{conf:?}");

        let pos: Vec<_> = (0..conf.len()).collect();
        let mut all_comps: Vec<bool> = vec![];
        for &i in &pos {
            let val = conf[i];
            //println!("i {i} = val {val}");
            let others: Vec<_> =
                pos.iter().filter(|&&j| j != i).map(|&j| conf[j]).collect();
            let pairs: Vec<_> = std::iter::repeat(val)
                .take(others.len())
                .zip(others)
                .collect();
            //dbg!(&pairs);
            let comps: Vec<_> =
                pairs.iter().map(|(x, y)| (x * 0.3) > *y as f64).collect();
            //dbg!(&comps);
            all_comps.push(comps.iter().all(|&v| v));
            println!(
                "{:8} {:0.04} {:5} => {:?}",
                families[i],
                val,
                comps.iter().all(|&v| v),
                comps
            );
        }
        //dbg!(&all_comps);
        let winners: Vec<String> = families
            .iter()
            .zip(all_comps)
            .filter_map(|(fam, win)| win.then(|| fam.to_string()))
            .collect();
        dbg!(&winners);

        if winners.len() == 1 {
            let winner = winners.first().unwrap();
            clear_winners
                .entry(winner.to_string())
                .and_modify(|v| *v += 1)
                .or_insert(1);
            println!("Clear Winner: {winner}");
        } else {
            println!("NO WINNER");
            let mut fam_comps: Vec<_> = conf.iter().zip(families).collect();
            fam_comps.sort_by(|a, b| {
                b.0.partial_cmp(a.0).unwrap().then_with(|| a.1.cmp(&b.1))
            });
            println!("{:?}", fam_comps);
            let (&top_conf, _) = fam_comps.first().unwrap();
            let threshold = top_conf / 5.;
            let winning_set: Vec<_> = fam_comps
                .iter()
                .filter_map(|(&c, f)| (c > threshold).then(|| f))
                .collect();
            println!("Winning Set: {winning_set:?}");
            for mut pair in winning_set.into_iter().permutations(2) {
                pair.sort();
                match pair[..] {
                    [&f1, &f2] => {
                        println!("f1 {f1} f2 {f2}");
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

    println!("Clear Winners");
    println!("{clear_winners:?}");

    println!("Winning Setts");
    println!("{winning_sets:?}");

    Ok(())
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
