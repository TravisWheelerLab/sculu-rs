use anyhow::{anyhow, bail, Result};
use clap::Parser;
use sculu::{Cli, ClusterArgs, Command, ComponentsArgs};
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
    let num_threads = args.num_threads.unwrap_or(num_cpus::get());

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
            sculu::build_components(args, num_threads)?;
            Ok(())
        }
        Some(Command::Cluster(args)) => {
            sculu::process_component(args, num_threads)?;
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
            let built_components = sculu::build_components(
                &ComponentsArgs {
                    alphabet: args.alphabet.clone(),
                    consensi: args.consensi.clone(),
                    instances: args.instances.clone(),
                    outdir: args.outdir.clone(),
                    config: args.config.clone(),
                },
                num_threads,
            )?;

            let mut merged: Vec<PathBuf> = vec![];
            for component in built_components.components {
                let res = sculu::process_component(
                    &ClusterArgs {
                        alphabet: args.alphabet.clone(),
                        consensi: built_components.consensi.clone(),
                        instances: built_components.instances_dir.clone(),
                        outdir: args.outdir.clone(),
                        config: args.config.clone(),
                        component,
                    },
                    num_threads,
                )?;

                merged.push(res);
            }

            let _ = sculu::concat_files(
                &built_components.consensi,
                &built_components.singletons,
                &merged,
                &args.outfile,
            )?;

            Ok(())
        }
        _ => unreachable!(),
    }
}
