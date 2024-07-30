use anyhow::Result;
use assert_cmd::Command;
use homedir::my_home;
use predicates::prelude::*;
use std::path::PathBuf;
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
    let args = [
        "--consensi",
        "tests/inputs/consensi.fa",
        "--instances",
        "tests/inputs/instances/AluY.fa",
        "tests/inputs/instances/AluYa5.fa",
        "tests/inputs/instances/AluYb8.fa",
        "tests/inputs/instances/AluYb9.fa",
        "tests/inputs/instances/AluYm1.fa",
        "--log",
        "debug",
        "--independence-threshold",
        ".8",
        "--confidence-margin",
        "3",
        "--threads",
        "8",
        "--alignment-matrix",
        MATRIX,
        "--rmblast-dir",
        &format!("{}/.local/bin", home.to_string_lossy().to_string()),
    ];
    run(&args)
}

// --------------------------------------------------
fn run(args: &[&str]) -> Result<()> {
    let output = Command::cargo_bin(PRG)?.args(args).output().unwrap();
    assert!(output.status.success());
    Ok(())
}
