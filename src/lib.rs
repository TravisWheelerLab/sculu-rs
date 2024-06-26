use anyhow::Result;
use clap::Parser;

/// SCULU subfamily clustering tool
#[derive(Debug, Parser)]
#[command(author, version, about)]
pub struct Args {
    /// FASTA file of subfamily consensus sequences
    #[arg(long, value_name = "CONSENSUS")]
    pub consensus: Option<String>,

    /// Directory of instance files for each subfamily
    #[arg(long, value_name = "INSTANCES")]
    pub instances: Option<String>,

    /// FASTA file of test sequences
    #[arg(long, value_name = "INSTANCES")]
    pub test_set: Option<String>,
}

// --------------------------------------------------
pub fn run(args: Args) -> Result<()> {
    dbg!(&args);
    Ok(())
}
