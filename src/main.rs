use clap::Parser;
use scrubby::utils::init_logger;
use scrubby::terminal::{App, Commands};


fn main() -> anyhow::Result<()> {
    
    let cli = App::parse();

    init_logger(cli.log_file);

    match &cli.command {
        Commands::Clean(args) => {
            args.validate_and_build()?.clean()?;
        },
        Commands::Classifer(args) => {
            args.validate_and_build()?.clean()?;
        },
        Commands::Alignment(args) => {
            args.validate_and_build()?.clean()?;
        },
        Commands::Download(args) => {
            if args.list {

            } else {
                args.validate_and_build()?.download_index()?;
            }
        },
        Commands::Difference(args) => {
            args.validate_and_build()?.compute()?;
        },
    }

    Ok(())
}