use anyhow::Result;
use env_logger::fmt::Color;
use env_logger::Builder;
use log::{Level, LevelFilter};
use metabuli::MetabuliSeqMode;
use scrub::ScrubbyTool;
use std::{collections::HashSet, io::Write, path::PathBuf};
use structopt::StructOpt;
use thiserror::Error;

mod align;
mod cli;
mod kraken;
mod scrub;
mod utils;
mod metabuli;
mod scrubby;

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
    cli.validate_kraken2()?;
    cli.validate_metabuli()?;

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
            reads,
            kraken_db,
            kraken_threads,
            kraken_taxa,
            kraken_taxa_direct,
            kraken_args,
            metabuli_db,
            metabuli_threads,
            metabuli_taxa,
            metabuli_taxa_direct,
            metabuli_seq_mode,
            metabuli_args,
            minimap2_index,
            minimap2_preset,
            minimap2_threads,
            minimap2_args,
            strobealign_index,
            strobealign_mode,
            strobealign_threads,
            strobealign_args,
            min_len,
            min_cov,
            min_mapq,
            output_format,
            compression_level,
        } => {

            let settings = scrub::Settings::new(
                [kraken_taxa.clone(), metabuli_taxa.clone()].concat(),
                [kraken_taxa_direct.clone(), metabuli_taxa_direct.clone()].concat(),
                min_len,
                min_cov,
                min_mapq,
                extract,
            );

            let mut scrubber =
                scrub::Scrubber::new(workdir, output_format, compression_level, settings, cli.force)?;

            scrubber.test_dependencies(
                !kraken_db.is_empty(), 
                !metabuli_db.is_empty(), 
                !minimap2_index.is_empty(), 
                !strobealign_index.is_empty()
            )?;

            let mut read_files = input;
            let mut reads_extract: HashSet<String> = HashSet::new();
            let mut scrub_index = 0;

            // Kraken2 taxonomic scrubbing
            match !kraken_db.is_empty() {
                false => log::info!("No databases specified: Kraken2"),
                true => {
                    
                    for db_path in kraken_db {
                        let db_name = get_reference_name(&db_path)?;

                        let (kraken_files, kraken_command) = scrubber.run_kraken(
                            &read_files,
                            &db_path,
                            &db_name,
                            &scrub_index,
                            &kraken_threads,
                            &kraken_args
                        )?;

                        let reads = scrubber.parse_kraken(
                            &kraken_files,
                            &kraken_taxa,
                            &kraken_taxa_direct,
                        )?;
                        scrubber.reads.add(&reads, ScrubbyTool::Kraken2, &db_name);

                        let (summary, files) = scrubber.deplete_to_workdir(
                            &read_files,
                            &reads,
                            ScrubbyTool::Kraken2,
                            &db_name,
                            db_path,
                            &scrub_index,
                            &extract,
                            &kraken_command
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


            // Metabuli taxonomic scrubbing
            match !metabuli_db.is_empty() {
                false => log::info!("No databases specified: Metabuli"),
                true => {
                    
                    for db_path in metabuli_db {
                        let db_name = get_reference_name(&db_path)?;

                        let (metabuli_files, metabuli_command) = scrubber.run_metabuli(
                            &read_files,
                            &db_path,
                            &db_name,
                            &scrub_index,
                            &metabuli_threads,
                            match metabuli_seq_mode {
                                None => None,
                                Some(ref seqmode) => Some(MetabuliSeqMode::from_arg(&seqmode))
                            },
                            &metabuli_args
                        )?;

                        let reads = scrubber.parse_metabuli(
                            &metabuli_files,
                            &metabuli_taxa,
                            &metabuli_taxa_direct,
                        )?;
                        scrubber.reads.add(&reads, ScrubbyTool::Metabuli, &db_name);

                        let (summary, files) = scrubber.deplete_to_workdir(
                            &read_files,
                            &reads,
                            ScrubbyTool::Metabuli,
                            &db_name,
                            db_path,
                            &scrub_index,
                            &extract,
                            &metabuli_command
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

                        let (alignment, minimap2_command) = scrubber.run_minimap2(
                            &read_files,
                            &index_path,
                            &index_name,
                            &scrub_index,
                            &minimap2_threads,
                            &minimap2_preset,
                            &minimap2_args
                        )?;

                        let reads = scrubber
                            .parse_alignment(&alignment, None, &min_len, &min_cov, &min_mapq)?;

                        scrubber.reads.add(&reads, ScrubbyTool::Minimap2, &index_name);

                        let (summary, files) = scrubber.deplete_to_workdir(
                            &read_files,
                            &reads,
                            ScrubbyTool::Minimap2,
                            &index_name,
                            index_path,
                            &scrub_index,
                            &extract,
                            &minimap2_command
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

                        let (alignment, strobealign_command) = scrubber.run_strobealign(
                            &read_files,
                            &index_path,
                            &index_name,
                            &scrub_index,
                            &strobealign_threads,
                            &strobealign_mode,
                            &strobealign_args
                        )?;

                        let reads = scrubber
                            .parse_alignment(&alignment, None, &min_len, &min_cov, &min_mapq)?;

                        scrubber.reads.add(&reads, ScrubbyTool::Strobealign, &index_name);

                        let (summary, files) = scrubber.deplete_to_workdir(
                            &read_files,
                            &reads,
                            ScrubbyTool::Strobealign,
                            &index_name,
                            index_path,
                            &scrub_index,
                            &extract,
                            &strobealign_command
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

                // This is inefficient - should be fixed but need for total counts...
                false => scrubber.write_depleted_pipeline_outputs(read_files, output)?,
            }

            scrubber.write_read_summary(reads)?;
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
            let krk_name = match kraken_name.clone() {
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
                scrub::Scrubber::new(workdir, output_format, compression_level, settings, cli.force)?;

            let reads = scrubber.parse_kraken(
                &Vec::from([kraken_report, kraken_reads.clone()]),
                &kraken_taxa,
                &kraken_taxa_direct,
            )?;
            scrubber.reads.add(&reads, ScrubbyTool::Kraken2, &match kraken_name.clone() { Some(name) => name.clone(), None => "kraken2".to_string()});

            let (summary, _) = scrubber.deplete_to_file(
                &input,
                &output,
                &reads,
                ScrubbyTool::Kraken2,
                &krk_name,
                kraken_reads,
                &0,
                &extract,
                &String::new()
            )?;

            scrubber.json.pipeline.push(summary.clone());
            scrubber.json.update(summary.total);
            scrubber.write_summary(json)?;
            scrubber.clean_up(false)?;
        }
        cli::Commands::ScrubMetabuli {
            input,
            output,
            workdir,
            extract,
            json,
            metabuli_report,
            metabuli_reads,
            metabuli_taxa,
            metabuli_taxa_direct,
            metabuli_name,
            output_format,
            compression_level,
        } => {
            let met_name = match metabuli_name.clone() {
                Some(name) => name,
                _ => get_reference_name(&metabuli_reads)?,
            };

            let settings = scrub::Settings::new(
                metabuli_taxa.clone(),
                metabuli_taxa_direct.clone(),
                0,
                0.,
                0,
                extract,
            );

            let mut scrubber =
                scrub::Scrubber::new(workdir, output_format, compression_level, settings, cli.force)?;
                        
            let reads = scrubber.parse_metabuli(
                &Vec::from([metabuli_report, metabuli_reads.clone()]),
                &metabuli_taxa,
                &metabuli_taxa_direct,
            )?;
            scrubber.reads.add(&reads, ScrubbyTool::Metabuli, &match metabuli_name { Some(name) => name, None => "metabuli".to_string()});

            let (summary, _) = scrubber.deplete_to_file(
                &input,
                &output,
                &reads,
                ScrubbyTool::Metabuli,
                &met_name,
                metabuli_reads,
                &0,
                &extract,
                &String::new()
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
            let aln_name = match alignment_name.clone() {
                Some(name) => name,
                _ => get_reference_name(&alignment)?,
            };

            let settings =
                scrub::Settings::new(Vec::new(), Vec::new(), min_len, min_cov, min_mapq, extract);

            let mut scrubber =
                scrub::Scrubber::new(workdir, output_format, compression_level, settings, cli.force)?;


            let reads = scrubber.parse_alignment(
                &alignment,
                alignment_format,
                &min_len,
                &min_cov,
                &min_mapq,
            )?;
            scrubber.reads.add(&reads, ScrubbyTool::Alignment, &match alignment_name { Some(name) => name, None => "alignment".to_string()});

            let (summary, _) = scrubber.deplete_to_file(
                &input, &output, &reads, ScrubbyTool::Alignment, &aln_name, alignment, &0, &extract, &String::new()
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
