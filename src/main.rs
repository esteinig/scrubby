use anyhow::Result;
use env_logger::fmt::Color;
use env_logger::Builder;
use log::{Level, LevelFilter};
use std::{collections::HashSet, io::Write, path::PathBuf};
use structopt::StructOpt;
use thiserror::Error;

mod align;
mod cli;
mod kraken;
mod scrub;
mod utils;

#[derive(Error, Debug)]
pub enum ScrubbyError {
    /// Represents a failure to convert a PathBuf converted OsStr into a String
    #[error("incorrect format of the reference path, are there non-UTF8 characters?")]
    InvalidReferencePath,
    /// Indicates a failure to obtain an absolute path
    #[error("reference name could not be obtained from {0}")]
    ReferenceNameExtraction(String),
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
            strobealign_index,
            strobealign_mode,
            strobealign_threads,
            min_len,
            min_cov,
            min_mapq,
            output_format,
            compression_level,
        } => {
            let settings = scrub::Settings::new(
                kraken_taxa.clone(),
                kraken_taxa_direct.clone(),
                min_len,
                min_cov,
                min_mapq,
                extract,
            );

            let mut scrubber =
                scrub::Scrubber::new(workdir, output_format, compression_level, settings)?;

            let mut read_files = input;
            let mut reads_extract: HashSet<String> = HashSet::new();
            let mut scrub_index = 0;

            // Kraken2 taxonomic scrubbing
            match !kraken_db.is_empty() {
                false => log::info!("No databases specified: Kraken2"),
                true => {
                    for db_path in kraken_db {
                        let db_name = get_reference_name(&db_path)?;

                        let kraken_files = scrubber.run_kraken(
                            &read_files,
                            &db_path,
                            &db_name,
                            &scrub_index,
                            &kraken_threads,
                        )?;

                        let reads = scrubber.parse_kraken(
                            &kraken_files,
                            &kraken_taxa,
                            &kraken_taxa_direct,
                        )?;

                        let (summary, files) = scrubber.deplete_to_workdir(
                            &read_files,
                            &reads,
                            Some("kraken2".to_string()),
                            &db_name,
                            db_path,
                            &scrub_index,
                            &extract,
                        )?;
                        scrubber.json.pipeline.push(summary);

                        match extract {
                            false => read_files = files, // update depleted intermediary files
                            true => reads_extract.extend(reads), // do not update intermediary files
                        }

                        scrub_index += 1
                    }
                }
            }

            // Minimap2 alignment scrubbing
            match !minimap2_index.is_empty() {
                false => log::info!("No indices specified: minimap2"),
                true => {
                    for index_path in minimap2_index {
                        let index_name = get_reference_name(&index_path)?;

                        let alignment = scrubber.run_minimap2(
                            &read_files,
                            &index_path,
                            &index_name,
                            &scrub_index,
                            &minimap2_threads,
                            &minimap2_preset,
                        )?;

                        let reads = scrubber
                            .parse_alignment(&alignment, None, &min_len, &min_cov, &min_mapq)?;

                        let (summary, files) = scrubber.deplete_to_workdir(
                            &read_files,
                            &reads,
                            Some("minimap2".to_string()),
                            &index_name,
                            index_path,
                            &scrub_index,
                            &extract,
                        )?;
                        scrubber.json.pipeline.push(summary);

                        match extract {
                            false => read_files = files, // update depleted intermediary files
                            true => reads_extract.extend(reads), // do not update intermediary files
                        }

                        scrub_index += 1;
                    }
                }
            }

            // Strobealign alignment scrubbing
            match !strobealign_index.is_empty() {
                false => log::info!("No indices specified: strobealign"),
                true => {
                    for index_path in strobealign_index {
                        let index_name = get_reference_name(&index_path)?;

                        let alignment = scrubber.run_strobealign(
                            &read_files,
                            &index_path,
                            &index_name,
                            &scrub_index,
                            &strobealign_threads,
                            &strobealign_mode,
                        )?;

                        let reads = scrubber
                            .parse_alignment(&alignment, None, &min_len, &min_cov, &min_mapq)?;

                        let (summary, files) = scrubber.deplete_to_workdir(
                            &read_files,
                            &reads,
                            Some("strobealign".to_string()),
                            &index_name,
                            index_path,
                            &scrub_index,
                            &extract,
                        )?;
                        scrubber.json.pipeline.push(summary);

                        match extract {
                            false => read_files = files, // update depleted intermediary files
                            true => reads_extract.extend(reads), // do not update intermediary files
                        }

                        scrub_index += 1;
                    }
                }
            }

            match extract {
                true => {
                    scrubber.write_extracted_pipeline_outputs(read_files, output, &reads_extract)?
                }
                false => scrubber.write_depleted_pipeline_outputs(read_files, output)?,
            }

            scrubber.write_summary(json)?;
            scrubber.clean_up(keep)?;
        }
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
            kraken_name,
            output_format,
            compression_level,
        } => {
            let krk_name = match kraken_name {
                Some(name) => name,
                _ => get_reference_name(&kraken_reads)?,
            };

            let settings = scrub::Settings::new(
                kraken_taxa.clone(),
                kraken_taxa_direct.clone(),
                0,
                0.,
                0,
                extract,
            );

            let mut scrubber =
                scrub::Scrubber::new(workdir, output_format, compression_level, settings)?;

            let reads = scrubber.parse_kraken(
                &Vec::from([kraken_report, kraken_reads.clone()]),
                &kraken_taxa,
                &kraken_taxa_direct,
            )?;

            let (summary, _) = scrubber.deplete_to_file(
                &input,
                &output,
                &reads,
                None,
                &krk_name,
                kraken_reads,
                &0,
                &extract,
            )?;

            scrubber.json.pipeline.push(summary.clone());
            scrubber.json.update(summary.total);
            scrubber.write_summary(json)?;
            scrubber.clean_up(false)?;
        }
        cli::Commands::ScrubAlignment {
            input,
            output,
            workdir,
            extract,
            json,
            alignment,
            alignment_name,
            alignment_format,
            min_len,
            min_cov,
            min_mapq,
            output_format,
            compression_level,
        } => {
            let aln_name = match alignment_name {
                Some(name) => name,
                _ => get_reference_name(&alignment)?,
            };

            let settings =
                scrub::Settings::new(Vec::new(), Vec::new(), min_len, min_cov, min_mapq, extract);

            let mut scrubber =
                scrub::Scrubber::new(workdir, output_format, compression_level, settings)?;

            let reads = scrubber.parse_alignment(
                &alignment,
                alignment_format,
                &min_len,
                &min_cov,
                &min_mapq,
            )?;

            let (summary, _) = scrubber.deplete_to_file(
                &input, &output, &reads, None, &aln_name, alignment, &0, &extract,
            )?;

            scrubber.json.pipeline.push(summary.clone());
            scrubber.json.update(summary.total);
            scrubber.write_summary(json)?;
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
    match db_path.file_stem() {
        Some(name) => Ok(name
            .to_os_string()
            .into_string()
            .map_err(|_| ScrubbyError::InvalidReferencePath)?),
        None => {
            return Err(ScrubbyError::ReferenceNameExtraction(format!(
                "{:?}",
                db_path
            )))
        }
    }
}

// Utility function to initialise the logger with color scheme
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
            orange_style
                .set_color(Color::Rgb(255, 102, 0))
                .set_bold(true);
            let mut apricot_style = buf.style();
            apricot_style
                .set_color(Color::Rgb(255, 195, 0))
                .set_bold(true);

            let msg = match record.level() {
                Level::Warn => (
                    orange_style.value(record.level()),
                    orange_style.value(record.args()),
                ),
                Level::Info => (
                    green_style.value(record.level()),
                    white_style.value(record.args()),
                ),
                Level::Debug => (
                    apricot_style.value(record.level()),
                    apricot_style.value(record.args()),
                ),
                Level::Error => (
                    red_style.value(record.level()),
                    red_style.value(record.args()),
                ),
                _ => (
                    white_style.value(record.level()),
                    white_style.value(record.args()),
                ),
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
