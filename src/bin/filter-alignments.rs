use anyhow::Result;
use clap::Parser;

/// Filter alignments
#[derive(Debug, Parser)]
#[command(author, version, about)]
pub struct Args {
    /// FASTA file of subfamily consensus sequences
    #[arg(value_name = "ALIGNMENTS")]
    pub alignments: String,

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
    dbg!(&args);

    //let mut output: Box<dyn Write> = match args.output.as_str() {
    //    "-" => Box::new(io::stdout()),
    //    out_name => Box::new(File::create(out_name)?),
    //};

    Ok(())
}
