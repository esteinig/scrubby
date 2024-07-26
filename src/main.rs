use clap::Parser;
use scrubby::prelude::*;
use scrubby::terminal::{App, Commands};


fn main() -> anyhow::Result<()> {
    
    let cli = App::parse();

    init_logger(cli.log_file);

    match &cli.command {
        Commands::Clean(args) => {
            args.validate_and_build()?.clean()?;
        },
    }

    Ok(())
}