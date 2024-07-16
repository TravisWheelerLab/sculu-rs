use anyhow::{anyhow, Result};
use clap::Parser;
use kseq::parse_reader;
use std::{
    collections::HashMap,
    fs::File,
    io::{self, BufRead, BufReader, Write},
};

/// Downsample
#[derive(Debug, Parser)]
#[command(author, version, about)]
pub struct Args {
    /// FASTA file
    #[arg(value_name = "FILE")]
    pub file: String,

    /// Number of sequences to keep
    #[arg(short, long, value_name = "NUM", default_value = "100")]
    pub number: usize,

    /// Output filename
    #[arg(short, long, value_name = "OUTPUT", default_value = "-")]
    pub output: String,
}

// --------------------------------------------------
fn main() {
    if let Err(e) = run(Args::parse()) {
        eprintln!("{e}");
        std::process::exit(1);
    }
}

// --------------------------------------------------
pub fn run(args: Args) -> Result<()> {
    let min_length = find_min_len(&args.file, &args.number)?;
    let mut output: Box<dyn Write> = match args.output.as_str() {
        "-" => Box::new(io::stdout()),
        out_name => Box::new(File::create(out_name)?),
    };
    let mut reader = parse_reader(open(&args.file)?)?;
    let mut taken = 0;

    while let Some(rec) = reader.iter_record()? {
        if args.number > 0 && taken == args.number {
            break;
        }

        let seq_len = rec.seq().len();
        if min_length > 0 && seq_len >= min_length {
            taken += 1;
            if rec.is_fasta() {
                writeln!(
                    output,
                    ">{}{}\n{}",
                    rec.head(),
                    rec.des(),
                    rec.seq()
                )?;
            } else {
                writeln!(
                    output,
                    "@{}{}\n{}\n{}\n{}",
                    rec.head(),
                    rec.des(),
                    rec.seq(),
                    if rec.sep().is_empty() { "+" } else { rec.sep() },
                    if rec.qual().is_empty() {
                        "-".repeat(rec.seq().len())
                    } else {
                        rec.qual().to_string()
                    },
                )?;
            }
        }
    }

    Ok(())
}

// --------------------------------------------------
fn find_min_len(filename: &str, number: &usize) -> Result<usize> {
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
        if &top_n >= number {
            break;
        }
    }

    Ok(min_len)
}

// --------------------------------------------------
fn open(filename: &str) -> Result<Box<dyn BufRead>> {
    match filename {
        "-" => Ok(Box::new(BufReader::new(io::stdin()))),
        _ => Ok(Box::new(BufReader::new(
            File::open(filename).map_err(|e| anyhow!("{}: {e}", filename))?,
        ))),
    }
}
