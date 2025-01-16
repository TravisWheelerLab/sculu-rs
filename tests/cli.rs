use anyhow::Result;
use assert_cmd::Command;
use homedir::my_home;
use predicates::prelude::*;
use std::path::PathBuf;
use tempfile::tempdir;

const PRG: &str = "sculu";
const MATRIX: &str = "tests/inputs/matrices/25p41g.matrix";

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
// TODO: getting "BLAST engine error: Error: Unknown error code 2"
// Run "make run" to verify instead
//#[test]
//fn run1() -> Result<()> {
//    let home = my_home()?.unwrap_or(PathBuf::from("/Users/kyclark"));
//    let args = vec![
//        "--consensi".to_string(),
//        "tests/inputs/consensi.fa".to_string(),
//        "--instances".to_string(),
//        "tests/inputs/instances/AluY.fa".to_string(),
//        "tests/inputs/instances/AluYa5.fa".to_string(),
//        "tests/inputs/instances/AluYb8.fa".to_string(),
//        "tests/inputs/instances/AluYb9.fa".to_string(),
//        "tests/inputs/instances/AluYm1.fa".to_string(),
//        "--independence-threshold".to_string(),
//        ".8".to_string(),
//        "--confidence-margin".to_string(),
//        "3".to_string(),
//        "--align-matrix".to_string(),
//        MATRIX.to_string(),
//        "--rmblast-dir".to_string(),
//        format!("{}/.local/bin", home.to_string_lossy().to_string()),
//    ];
//    run(args)
//}
//
//// --------------------------------------------------
//fn run(mut args: Vec<String>) -> Result<()> {
//    let outdir = tempdir()?;
//    let outfile = outdir.path().join("final.fa");
//
//    args.extend_from_slice(&[
//        "--outfile".to_string(),
//        outfile.to_string_lossy().to_string(),
//        "--outdir".to_string(),
//        outdir.path().to_string_lossy().to_string(),
//    ]);
//
//    dbg!(Command::cargo_bin(PRG)?.args(args.clone()));
//
//    let output = Command::cargo_bin(PRG)?.args(args).output().unwrap();
//    dbg!(&output);
//    assert!(output.status.success());
//
//    // Ensure the output directory was created
//    assert!(outdir.path().is_dir());
//
//    // Ensure the final consensi file was created
//    assert!(outfile.exists());
//
//    Ok(())
//}
