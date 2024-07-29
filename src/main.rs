use clap::Parser;
use scrubby::utils::init_logger;
use scrubby::terminal::{App, Commands};


fn main() -> anyhow::Result<()> {
    
    let cli = App::parse();

    init_logger(cli.log_file);

    match &cli.command {
        Commands::Reads(args) => {
            args.validate_and_build()?.clean()?;
        },
        Commands::Classifer(args) => {
            args.validate_and_build()?.clean()?;
        },
        Commands::Alignment(args) => {
            args.validate_and_build()?.clean()?;
        },
        Commands::Download(args) => {
            let dl = args.validate_and_build()?;
            if args.list { dl.list() } else { dl.download_index()? };
        },
        Commands::Diff(args) => {
            args.validate_and_build()?.compute()?;
        },
    }

    Ok(())
}