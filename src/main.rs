use clap::Parser;
use sculu::Args;

// --------------------------------------------------
fn main() {
    if let Err(e) = sculu::run(Args::parse()) {
        eprintln!("Error: {e}");
        std::process::exit(1);
    }
}
