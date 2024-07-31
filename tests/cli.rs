use anyhow::Result;
use assert_cmd::Command;
use homedir::my_home;
use predicates::prelude::*;
use std::ffi::OsStr;
use std::path::PathBuf;
use tempfile::tempdir;
use walkdir::WalkDir;
//use pretty_assertions::assert_eq;

const PRG: &str = "sculu";
const MATRIX: &str = "tests//inputs/matrices/25p41g.matrix";

// --------------------------------------------------
#[test]
fn usage() -> Result<()> {
    for flag in &["-h", "--help"] {
        Command::cargo_bin(PRG)?
            .arg(flag)
            .assert()
            .stdout(predicate::str::contains("Usage"));
    }
    Ok(())
}

// --------------------------------------------------
#[test]
fn run1() -> Result<()> {
    let home = my_home()?.unwrap_or(PathBuf::from("/Users/kyclark"));
    let args = vec![
        "--consensi".to_string(),
        "tests/inputs/consensi.fa".to_string(),
        "--instances".to_string(),
        "tests/inputs/instances/AluY.fa".to_string(),
        "tests/inputs/instances/AluYa5.fa".to_string(),
        "tests/inputs/instances/AluYb8.fa".to_string(),
        "tests/inputs/instances/AluYb9.fa".to_string(),
        "tests/inputs/instances/AluYm1.fa".to_string(),
        "--log".to_string(),
        "debug".to_string(),
        "--independence-threshold".to_string(),
        ".8".to_string(),
        "--confidence-margin".to_string(),
        "3".to_string(),
        "--threads".to_string(),
        "8".to_string(),
        "--alignment-matrix".to_string(),
        MATRIX.to_string(),
        "--rmblast-dir".to_string(),
        format!("{}/.local/bin", home.to_string_lossy().to_string()),
    ];
    run(args)
}

// --------------------------------------------------
fn run(mut args: Vec<String>) -> Result<()> {
    let outdir = tempdir()?;
    args.extend_from_slice(&[
        "--outdir".to_string(),
        outdir.path().to_string_lossy().to_string(),
    ]);
    let output = Command::cargo_bin(PRG)?.args(args).output().unwrap();
    assert!(output.status.success());

    // Ensure the output directory was created
    assert!(outdir.path().is_dir());

    // Ensure the final consensi file was created
    let wanted = OsStr::new("final-consensi.fa");
    let files = WalkDir::new(outdir.path())
        .into_iter()
        .filter_map(Result::ok)
        .filter(|entry| {
            entry.file_type().is_file() && entry.file_name() == wanted
        })
        .collect::<Vec<_>>();

    assert_eq!(files.len(), 1);

    Ok(())
}
