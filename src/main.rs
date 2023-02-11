use anyhow::Result;
use thiserror::Error;
use structopt::StructOpt;
use std::path::PathBuf;
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
            output_format,
            compression_level,
        } => {
            
            let scrubber = scrub::Scrubber::new(workdir)?;

            log::info!("Welcome to Scrubby, your trusty read scrubber!");
            scrubber.run_kraken(&input, kraken_db, kraken_threads);

        }
    }

    log::info!("Scrub scrub, scrubbity-scrub! Your sequence data, only cleaner!");

    Ok(())
}