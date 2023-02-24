use anyhow::Result;
use thiserror::Error;
use structopt::StructOpt;
use env_logger::Builder;
use log::{LevelFilter, Level};
use env_logger::fmt::Color;
use std::{io::Write, path::PathBuf};

mod cli;
mod scrub;
mod kraken;
mod align;
mod utils;



#[derive(Error, Debug)]
pub enum ScrubbyError {
    /// Represents a failure to convert a PathBuf converted OsStr into a String
    #[error("incorrect format of the reference path, are there non-UTF8 characters?")]
    InvalidReferencePath,
    /// Indicates a failure to obtain an absolute path
    #[error("reference name could not be obtained from {0}")]
    ReferenceNameExtraction(String)
}

fn main() -> Result<()> {
    let cli = cli::Cli::from_args();

    // Additional command-line client checks
    cli.validate_input_output_combination()?;

    // Initiate color-schemed logger
    init_logger()?;

    log::info!("=============================================");
    log::info!("Welcome to Scrubby! You name it, we clean it.");
    log::info!("=============================================");

    match cli.commands {
        cli::Commands::ScrubReads {
            input,
            output,
            workdir,
            extract,
            keep,
            json,
            kraken_db,
            kraken_threads,
            kraken_taxa,
            kraken_taxa_direct,
            minimap2_index,
            minimap2_preset,
            minimap2_threads,
            min_len,
            min_cov,
            min_mapq,
            output_format,
            compression_level,
        } => {
            
            let mut scrubber = scrub::Scrubber::new(workdir, output_format, compression_level)?;
            
            let mut read_files = input;
            let mut scrub_index = 0;

            // Kraken2 taxonomic scrubbing
            match kraken_db.len() > 0 {
                false => log::info!("No databases specified: Kraken2"),
                true => {
                    for db_path in kraken_db.into_iter() {
                        let db_name = get_reference_name(&db_path)?;
                        let kraken_files = scrubber.run_kraken(&read_files, &db_path, &db_name, &scrub_index, &kraken_threads)?;
                        let (depletion_summary, files) = scrubber.deplete_kraken(
                            &read_files,
                            None,
                            &db_name,
                            &scrub_index,
                            &false,
                            &kraken_files,
                            &kraken_taxa,
                            &kraken_taxa_direct
                        )?;
                        read_files = files;
                        scrubber.json.summary.push(depletion_summary);
                        scrub_index += 1
                    }
                }
            }

            // Minimap2 alignment scrubbing
            match minimap2_index.len() > 0 {
                false => log::info!("No indices specified: minimap2"),
                true => {
                    for index_path in minimap2_index.into_iter() {
                        let index_name = get_reference_name(&index_path)?;
                        let alignment = scrubber.run_minimap2(&read_files, &index_path, &index_name, &scrub_index, &minimap2_threads, &minimap2_preset)?;
                        let (depletion_summary, files) = scrubber.deplete_alignment(
                            &read_files, 
                            None,
                            &alignment, 
                            None,
                            &index_name,
                            &scrub_index,
                            &false,
                            &min_len,
                            &min_cov,
                            &min_mapq
                        )?;
                        read_files = files;
                        scrubber.json.summary.push(depletion_summary);
                        scrub_index += 1;
                    }
                }
            }
            

            // Iterating again over the depleted record files to produce the user-specified outputs
            // because we want to ensure they are properly compressed (intermediate files are not)
            scrubber.write_outputs(read_files, output)?;

            // Summary output depending on the value of --json: None, file path or "-"
            scrubber.write_summary(json)?;
            
            // If we do not want to keep the intermediary files in `workdir` delete the directory
            scrubber.clean_up(keep)?;

        },
        cli::Commands::ScrubKraken { 
            input,
            output,
            workdir,
            extract,
            json,
            kraken_report,
            kraken_reads,
            kraken_taxa,
            kraken_taxa_direct,
            output_format,
            compression_level
         } => {

            let mut scrubber = scrub::Scrubber::new(workdir, output_format, compression_level)?;

            let kraken_name = get_reference_name(&kraken_reads)?;
            let (depletion_summary, _) = scrubber.deplete_kraken(
                &input,
                Some(output),
                &kraken_name,
                &0, 
                &false,
                &Vec::from([kraken_report, kraken_reads]),
                &kraken_taxa,
                &kraken_taxa_direct
            )?;
            
            scrubber.json.summary.push(depletion_summary);

            // Summary output depending on the value of --json: None, file path or "-"
            scrubber.write_summary(json)?;
            
            // A bit inefficient but we can fix the workdir creation later, 
            // for now always delete since we do not actually use it
            scrubber.clean_up(false)?;

         }
        cli::Commands::ScrubAlignment { 
            input,
            output,
            workdir,
            extract,
            json,
            alignment,
            alignment_format,
            min_len,
            min_cov,
            min_mapq,
            output_format,
            compression_level
        } => {

            let mut scrubber = scrub::Scrubber::new(workdir, output_format, compression_level)?;
            let alignment_name = get_reference_name(&alignment)?;
            let (depletion_summary, _) =  scrubber.deplete_alignment(
                &input, 
                Some(output),
                &alignment, 
                alignment_format,
                &alignment_name,
                &0,
                &false,
                &min_len,
                &min_cov,
                &min_mapq
            )?;
            
            scrubber.json.summary.push(depletion_summary);
        
            // Summary output depending on the value of --json: None, file path or "-"
            scrubber.write_summary(json)?;
            
            // A bit inefficient but we can fix the workdir creation later, 
            // for now always delete since we do not actually use it
            scrubber.clean_up(false)?;

        }
    }
    log::info!("==============================================================");
    log::info!("Thank you for using Scrubby! Your sequence data, only cleaner.");
    log::info!("==============================================================");

    Ok(())
}

// Utility function to extract the database name as valid UTF-8
fn get_reference_name(db_path: &PathBuf) -> Result<String, ScrubbyError> {
    match db_path.file_stem(){
        Some(name) => Ok(name.to_os_string().into_string().map_err(|_| ScrubbyError::InvalidReferencePath)?),
        None => return Err(ScrubbyError::ReferenceNameExtraction(format!("{:?}", db_path)))
    }
}

//Utility function to initialise the logger with color scheme
fn init_logger() -> Result<()> {
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
    Ok(())
}