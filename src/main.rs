use clap::Parser;
use sculu::{self, Args};

// --------------------------------------------------
fn main() {
    if let Err(e) = sculu::run(Args::parse()) {
        eprintln!("{e}");
        std::process::exit(1);
    }
}
