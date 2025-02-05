use anyhow::Result;
use clap::Parser;
use sculu::{self, Cli, Command};
use std::{fs::{self, File}, io::BufWriter};
use tempfile::tempdir;

// --------------------------------------------------
fn main() {
    if let Err(e) = run(Cli::parse()) {
        eprintln!("Error: {e}");
        std::process::exit(1);
    }
}

// --------------------------------------------------
fn run(cli_args: Cli) -> Result<()> {
    // Use a tempdir if no named outdir
    let outdir = cli_args
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

    match &cli_args.command {
        Some(Command::Cluster(cluster_args)) => {
            sculu::cluster(&outdir, &cli_args, cluster_args.clone())?;
            Ok(())
        }
        Some(Command::Merge(merge_args)) => {
            sculu::merge(&outdir, &cli_args, merge_args.clone())?;
            Ok(())
        }
        _ => unreachable!(),
    }
}
