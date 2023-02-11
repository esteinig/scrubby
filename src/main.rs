use anyhow::Result;
use structopt::StructOpt;
use chrono::Local;
use env_logger::Builder;
use log::LevelFilter;
use std::io::Write;
mod scrub;
mod cli;

fn main() -> Result<()> {
    let cli = cli::Cli::from_args();

    // Command specific checks - scrubbing

    cli.validate_input_output_combination()?;

    Builder::new()
        .format(|buf, record| {
            writeln!(buf,
                "{} [{}] - {}",
                Local::now().format("%Y-%m-%dT%H:%M:%S"),
                record.level(),
                record.args()
            )
        })
        .filter(None, LevelFilter::Info)
        .init();


    match cli.commands {
        cli::Commands::Scrub {
            input,
            output,
            workdir,
            kraken_db,
            kraken_threads,
            kraken_taxa,
            kraken_taxa_direct,
            output_format,
            compression_level,
        } => {
            
            let scrubber = scrub::Scrubber::new(workdir)?;

            log::info!("Welcome to Scrubby! You name it, we clean it.");
            for db in kraken_db{
                scrubber.run_kraken(&input, db, kraken_threads);
                scrubber.deplete_kraken(&input, &kraken_taxa, &kraken_taxa_direct)?;
            }
            
        }
    }

    log::info!("Thank you for using Scrubby.");
    log::info!("Your sequence data, only cleaner.");

    Ok(())
}