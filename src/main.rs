use anyhow::Result;
use thiserror::Error;
use structopt::StructOpt;
use env_logger::Builder;
use log::{LevelFilter, Level};
use env_logger::fmt::Color;
use std::{io::Write, path::PathBuf};

mod scrub;
mod kraken;
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

    // Command line application specific checks
    cli.validate_input_output_combination()?;

    Builder::new()
        .format(|buf, record| {
            let timestamp = buf.timestamp();

            let mut red_style = buf.style();
            red_style.set_color(Color::Red).set_bold(true);
            let mut green_style = buf.style();
            green_style.set_color(Color::Green).set_bold(true);
            let mut white_style = buf.style();
            white_style.set_color(Color::White).set_bold(false);
            let mut orange_style = buf.style();
            orange_style.set_color(Color::Rgb(255, 102, 0)).set_bold(true);
            let mut apricot_style = buf.style();
            apricot_style.set_color(Color::Rgb(255, 195, 0)).set_bold(true);

            let msg = match record.level(){
                Level::Warn => (orange_style.value(record.level()), orange_style.value(record.args())),
                Level::Info => (green_style.value(record.level()), white_style.value(record.args())),
                Level::Debug => (apricot_style.value(record.level()), apricot_style.value(record.args())),
                Level::Error => (red_style.value(record.level()), red_style.value(record.args())),
                _ => (white_style.value(record.level()), white_style.value(record.args()))
            };

            writeln!(
                buf,
                "{} [{}] - {}",
                white_style.value(timestamp),
                msg.0,
                msg.1
            )
        })
        .filter(None, LevelFilter::Info)
        .init();


    match cli.commands {
        cli::Commands::ScrubReads {
            input,
            output,
            workdir,
            extract,
            keep,
            kraken_db,
            kraken_threads,
            kraken_taxa,
            kraken_taxa_direct,
            output_format,
            compression_level,
        } => {
            
            let scrubber = scrub::Scrubber::new(workdir, output_format, compression_level)?;
            
            log::info!("=============================================");
            log::info!("Welcome to Scrubby! You name it, we clean it.");
            log::info!("=============================================");

            let mut read_files = input;
            for (db_index, db_path) in kraken_db.into_iter().enumerate() {
                // Extract the database name from path and run Kraken2
                let db_name = get_db_name(&db_path)?;
                let kraken_files = scrubber.run_kraken(&read_files, &db_path, &db_name, &db_index, &kraken_threads)?;
                // These are either depleted or extracted reads - for depleted reads, we use the depleted reads as input for the next iteration
                // but for extracted reads, we do not want to re-extract reads with another database [https://github.com/esteinig/scrubby/issues/2]
                read_files = scrubber.deplete_kraken(&read_files, &db_name, &db_index, &false, &kraken_files, &kraken_taxa, &kraken_taxa_direct)?;
            }
            // Iterating again over the depleted record files to produce the user-specified outputs
            // because we want to ensure they are properly compressed (intermediate files are not)
            scrubber.write_outputs(read_files, output)?;
            // If we do not want to keep the intermediary files in `workdir` delete the directory
            scrubber.clean_up(keep)?;


            
        }
    }
    log::info!("==============================================================");
    log::info!("Thank you for using Scrubby! Your sequence data, only cleaner.");
    log::info!("==============================================================");

    Ok(())
}

// Utility function to extract the database name as valid UTF-8
fn get_db_name(db_path: &PathBuf) -> Result<String, ScrubbyError> {
    match db_path.file_name(){
        Some(name) => Ok(name.to_os_string().into_string().map_err(|_| ScrubbyError::InvalidDatabasePath)?),
        None => return Err(ScrubbyError::DatabaseNameExtraction(format!("{:?}", db_path)))
    }
}