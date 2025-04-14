use anyhow::{anyhow, bail, Result};
use clap::Parser;
use sculu::{Cli, Command};
use std::{
    fs::{self, File},
    io::BufWriter,
    path::PathBuf,
};

// --------------------------------------------------
fn main() {
    if let Err(e) = run(Cli::parse()) {
        eprintln!("Error: {e}");
        std::process::exit(1);
    }
}

// --------------------------------------------------
fn run(args: Cli) -> Result<()> {
    let log_file = args.logfile.unwrap_or(PathBuf::from("debug.log"));

    if let Some(dir) = log_file.parent() {
        if !dir.exists() {
            fs::create_dir_all(dir)?;
        }
    }

    env_logger::Builder::new()
        .filter_level(log::LevelFilter::Debug)
        .target(match log_file.to_str() {
            Some("-") => env_logger::Target::Stdout,
            Some(path) => env_logger::Target::Pipe(Box::new(BufWriter::new(
                File::create(&log_file).map_err(|e| anyhow!("{path}: {e}"))?,
            ))),
            _ => bail!(r#"Cannot parse --log-file "{}""#, log_file.display()),
        })
        .init();

    match &args.command {
        Some(Command::Config(args)) => {
            sculu::write_config(args)?;
            Ok(())
        }
        Some(Command::Components(args)) => {
            sculu::build_components(args)?;
            Ok(())
        }
        Some(Command::Cluster(args)) => {
            sculu::process_component(args)?;
            Ok(())
        }
        Some(Command::Concat(args)) => {
            sculu::concat_files(
                &args.consensi,
                &args.singletons,
                &args.components,
                &args.outfile,
            )?;
            Ok(())
        }
        Some(Command::Run(args)) => {
            let components = sculu::build_components(args)?;
            let mut merged: Vec<PathBuf> = vec![];
            for component in components.components {
                let mut process_args = args.clone();
                process_args.component = Some(component);
                let res = sculu::process_component(&process_args)?;
                merged.push(res);
            }

            let _ = sculu::concat_files(
                &args.consensi,
                &components.singletons,
                &merged,
                &args.outfile,
            )?;

            Ok(())
        }
        _ => unreachable!(),
    }
}
