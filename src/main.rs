use clap::Parser;
use sculu::{Args, run};

// --------------------------------------------------
fn main() {
    if let Err(e) = run(Args::parse()) {
        eprintln!("Error: {e}");
        std::process::exit(1);
    }
}
