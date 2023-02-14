use anyhow::Result;
use thiserror::Error;
use structopt::StructOpt;
use chrono::Local;
use env_logger::Builder;
use log::LevelFilter;
use std::{io::Write, ffi::OsString, path::PathBuf};

mod scrub;
mod cli;
mod utils;



#[derive(Error, Debug)]
pub enum ScrubbyError {
    /// Represents a failure to convert a PathBuf converted OsStr into a String
    #[error("incorrect format of the input database path, are there non-UTF8 characters?")]
    InvalidDatabasePath,
    /// Indicates a failure to obtain an absolute path
    #[error("database name could not be obtained from {0}")]
    DatabaseNameExtraction(String),
}

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
            extract,
            kraken_db,
            kraken_threads,
            kraken_taxa,
            kraken_taxa_direct,
            output_format,
            compression_level,
        } => {
            
            let scrubber = scrub::Scrubber::new(workdir, output_format, compression_level)?;

            log::info!("Welcome to Scrubby! You name it, we clean it.");

            let mut scrubbed_reads = input;
            for (db_index, db_path) in kraken_db.into_iter().enumerate() {

                let db_name = get_db_name(&db_path)?;

                let kraken_files = scrubber.run_kraken(&scrubbed_reads, &db_path, &db_name, &db_index, &kraken_threads)?;
                // These are either depleted or extracted reads - for depleted reads, we use the depleted reads as input for the next iteration
                // but for extracted reads, we do not want to re-extract reads with another database [https://github.com/esteinig/scrubby/issues/2]
                scrubbed_reads = scrubber.deplete_kraken(&scrubbed_reads, &db_name, &db_index, &false, &kraken_files, &kraken_taxa, &kraken_taxa_direct)?;
            }
            
        }
    }

    log::info!("Thank you for using Scrubby! Your sequence data, only cleaner.");

    Ok(())
}

fn get_db_name(db_path: &PathBuf) -> Result<String, ScrubbyError> {
    match db_path.file_name(){
        // Kraken2 database name conversion - database path must be UTF8-convertable for manipulation of file paths
        Some(name) => Ok(name.to_os_string().into_string().map_err(|_| ScrubbyError::InvalidDatabasePath)?),
        None => return Err(ScrubbyError::DatabaseNameExtraction(format!("{:?}", db_path)))
    }
}