use chrono::Local;
use std::fs::File;
use anyhow::Result;
use serde::Serialize;
use thiserror::Error;
use std::path::PathBuf;
use std::str::from_utf8;
use std::process::Command;
use std::fs::{create_dir_all, remove_dir_all};
use std::collections::HashSet; 
use needletail::{parse_fastx_file, FastxReader};
use std::io::BufWriter;

use crate::utils::CompressionExt;


#[derive(Error, Debug)]
pub enum ScrubberError {
    /// Represents a failure to execute Kraken2
    #[error("failed to run Kraken2 - is it installed?")]
    KrakenExecutionError,
    /// Represents a failure to run Kraken2
    #[error("failed to run Kraken2")]
    KrakenClassificationError,
    /// Represents a failure to count a taxonomic parent during report parsing
    #[error("failed to provide a parent taxon while parsing report from Kraken2")]
    KrakenReportTaxonParent,
    /// Represents a failure to convert the read field from string to numeric field in the report file
    #[error("failed to convert the read field in the report from Kraken2")]
    KrakenReportReadFieldConversion,
    /// Represents a failure to convert the read field from string to numeric field in the report file
    #[error("failed to convert the direct read field in the report from Kraken2")]
    KrakenReportDirectReadFieldConversion,
    /// Represents a failure in the correct length of the input file vector
    #[error("input file error - incorrect number of input files")]
    FileNumberError,
    /// Represents a failure to convert a PathBuf converted OsStr into a String
    #[error("input file error - incorrect format of the input file path: are there non-standard characters?")]
    InvalidFilePathConversion,
    /// Represents all other cases of `std::io::Error`.
    #[error(transparent)]
    IOError(#[from] std::io::Error),
    /// Indicates that the working directory already exists
    #[error("working directory exists at {0}")]
    WorkdirExists(String),
    /// Indicates a failure to create the working directory
    #[error("working directory path could not be created from {0}")]
    WorkdirCreate(String),
    /// Indicates a failure to obtain an absolute path
    #[error("absolute path could not be obtained from {0}")]
    AbsolutePath(String),
    /// Indicates failure to parse file with Needletail
    #[error("failed file input/output")]
    FastxRecordIO(#[source] needletail::errors::ParseError),
    /// Indicates failure to obntain compression writer with Niffler
    #[error("failed to get compression writer")]
    DepletionCompressionWriter(#[source] niffler::Error),
    /// Indicates failure to parse record identifier
    #[error("failed to parse record identifier")]
    DepletionRecordIdentifier(#[source] std::str::Utf8Error),
}

pub struct Scrubber {
    workdir: PathBuf,
    output_format: Option<niffler::compression::Format>,
    compression_level: niffler::compression::Level 
}

impl Scrubber {
    pub fn new(workdir: Option<PathBuf>, output_format: Option<niffler::compression::Format>, compression_level: niffler::compression::Level ) -> Result<Self, ScrubberError> {
        let _workdir = check_or_create_workdir(workdir)?;
        Ok(Self { workdir: _workdir, output_format, compression_level })
    }
    ///
    pub fn run_kraken(
        &self,
        input: &Vec<PathBuf>,
        db_path: &PathBuf,
        db_name: &String,
        db_index: &usize,
        threads: &u32,
    ) -> Result<Vec<PathBuf>, ScrubberError>{
        
        // Safely build the arguments for Kraken2
        let kraken_args = crate::kraken::get_kraken_command(input, db_path, db_name, db_index, threads)?;

        log::info!("Executing taxonomic classification with Kraken2 ({})", db_name);
        log::debug!("Executing Kraken2 command: {}", &kraken_args.join(" "));

        // Run the Kraken command
        let output = Command::new("kraken2")
            .args(kraken_args)
            .current_dir(&self.workdir)
            .output()
            .map_err(|_| ScrubberError::KrakenExecutionError)?;

        // Ensure command ran successfully
        if output.status.success() {
            log::info!("Completed taxonomic classification with Kraken2 ({})", db_name)
        } else {
            log::error!("Failed to run taxonomic classification with Kraken2 ({})", db_name);
            let err_msg = crate::kraken::get_kraken_err_msg(output)?;
            log::error!("Error from {}", err_msg);
            return Err(ScrubberError::KrakenClassificationError)
        }

        let kraken_report = self.workdir.join(format!("{}-{}.report", db_index, db_name));
        let kraken_reads = self.workdir.join(format!("{}-{}.kraken", db_index, db_name));
        
        Ok(Vec::from([kraken_report, kraken_reads]))
    }
    ///
    pub fn deplete_kraken(
        &self,
        input: &Vec<PathBuf>,
        db_name: &String,
        db_index: &usize,
        extract: &bool,
        kraken_files: &Vec<PathBuf>,
        kraken_taxa: &Vec<String>,
        kraken_taxa_direct: &Vec<String>
    ) -> Result<Vec<PathBuf>, ScrubberError>{

        log::info!("Parsing taxonomic classification report...");
        
        let taxids = crate::kraken::get_taxids_from_report(
            kraken_files[0].clone(),
            kraken_taxa,
            kraken_taxa_direct
        )?;

        let reads = crate::kraken::get_taxid_reads(
            taxids, kraken_files[1].clone()
        )?;
        
        // Temporary files for sequential depletion/extraction in workdir
        let output = match input.len() {
                2 => Vec::from([self.workdir.join(format!("{}-{}_1.fq", db_index, db_name)),  self.workdir.join(format!("{}-{}_2.fq", db_index, db_name))]),
                1 => Vec::from([self.workdir.join(format!("{}-{}.fq", db_index, db_name))]),
                _ => return Err(ScrubberError::FileNumberError)
            };

        let mut read_summary = ReadCountSummary::new();

        // Enumerated loop is ok since we have checked matching input/output 
        // vector length in command-line interface and match the file number
        for (i, _) in input.iter().enumerate() {
            // Initiate the depletion operator and deplete/extract the reads identifiers parsed from Kraken
            let depletor = ReadDepletion::new(self.output_format, self.compression_level)?;
            let read_counts = depletor.deplete(&reads, &input[i], &output[i], extract)?;

            read_summary.reads.push(read_counts);
        };

        Ok(output)
    }
    /// 
    pub fn write_outputs(&self, input_files: Vec<PathBuf>, output_files: Vec<PathBuf>) -> Result<(), ScrubberError> {

        for (i, _) in input_files.iter().enumerate() { // input output are iensured to have same length through command-line interface
            let (mut reader, mut writer) = get_reader_writer(&input_files[i], &output_files[i], self.compression_level, self.output_format)?;

            log::info!("Writing reads to output file: {:?}", output_files[i]);

            while let Some(record) = reader.next() {
                let rec = record.map_err(|err| ScrubberError::FastxRecordIO(err))?;
                rec.write(&mut writer, None).map_err(|err| ScrubberError::FastxRecordIO(err))?;
            }
        }   
        Ok(())
    }
    pub fn clean_up(&self, keep: bool) -> Result<(), ScrubberError> {
        match keep {
            true => {
                log::info!("Keeping working directory and intermediary files");
            },
            false => {
                log::info!("Deleting working directory and intermediary files");
                remove_dir_all(&self.workdir)?;
            }
        }
        Ok(())
    }
}



/// Checks if work directory exists and otherwise creates it
///  
/// # Errors
/// A [`ScrubberError::WorkdirExists`](#scrubbererror) is returned if the directory already exists
/// A [`ScrubberError::WorkdirCreate`](#scrubbererror) is returned if the directory cannot be created
/// A [`ScrubberError::AbsolutePath`](#scrubbererror) is returned if the directory path cannot be canonicalized
pub fn check_or_create_workdir(workdir: Option<PathBuf>) -> Result<PathBuf, ScrubberError> {
    let _workdir = match workdir {
        Some(path) => path,
        None => PathBuf::from(format!("Scrubby-{}", Local::now().format("%Y%m%dT%H%M%S")))
    };
    let workdir_msg = format!("{:?}", _workdir);
    if !&_workdir.exists(){
        create_dir_all(&_workdir).map_err(|_| ScrubberError::WorkdirCreate(workdir_msg))?;
        let abs_workdir = _workdir.canonicalize().map_err(|_| ScrubberError::AbsolutePath(format!("{:?}", _workdir)))?;
        Ok(abs_workdir)
    } else {
        Err(ScrubberError::WorkdirExists(workdir_msg))
    }
}


/*
=========================
Read depletion/extraction
=========================
*/

// A summary struct to hold counts
// for the read depletion/extraction
#[derive(Serialize)]
pub struct ReadCounts {
    pub total: u64,
    pub depleted: u64,
    pub extracted: u64,
    pub retained: u64,
    pub input_file: PathBuf,
    pub output_file: PathBuf,
}

// A summary struct to hold settings
// for depletion/extraction
#[derive(Serialize)]
pub struct DepletionSettings {
}

// Struct to hold the read depletion for
// each file to be output to JSON
#[derive(Serialize)]
pub struct ReadCountSummary {
    pub reads: Vec<ReadCounts>,
    pub settings: DepletionSettings,
}

impl ReadCountSummary {
    pub fn new() -> Self {
        Self {
            reads: Vec::new(),
            settings: DepletionSettings {},
        }
    }
}

/// Read depletion struct
/// 
/// This struct can be initialised with `niffler ::compression::Format` which
/// optionally specifies the output format and a `niffler::compression::Level`
/// which specified the output compression level.
pub struct ReadDepletion {
    output_format: Option<niffler::compression::Format>,
    compression_level: niffler::compression::Level
}

impl ReadDepletion {
    pub fn new(
        output_format: Option<niffler::compression::Format>,
        compression_level: niffler::compression::Level
    ) -> Result<Self, ScrubberError> {
                
        Ok(Self { output_format, compression_level })
    }
    /// Deplete reads from an input read file
    /// 
    /// This method depletes (`extract = false`) or extracts (`extract = true`) an input
    /// read file which may be compressed. It checks if the read identifier is contained
    /// within the `HashSet`of read identifiers provided. It also counts the total,
    /// depleted or extracted, and the retained reads and returns a `ReadCounts` object
    /// which can be added to the `ReadCountSummary` to be output to JSON.
    /// 
    /// /// # Errors
    /// A [`ScrubberError::DepletionFastxParser`](#scrubbererror) is returned if the input file cannot be read, or if the output file cannot be written to
    /// A [`ScrubberError::DepletionCompressionWriter`](#scrubbererror) is returned if the compression writer cannot be obtained
    /// A [`ScrubberError::DepletionRecordIdentifier`](#scrubbererror) if the read identifier cannot be converted to valid UTF8
    /// 
    pub fn deplete(
        &self, 
        reads: &HashSet<String>,
        input: &PathBuf,
        output: &PathBuf,
        extract: &bool,
    ) -> Result<ReadCounts, ScrubberError> {

        // Input output of read files includes compression detection
        let (mut reader, mut writer) = get_reader_writer(input, output, self.compression_level, self.output_format)?;

        let mut read_counts = ReadCounts {
            total: 0,
            depleted: 0,
            extracted: 0,
            retained: 0,
            input_file: input.to_path_buf(),
            output_file: output.to_path_buf(),
        };
        
        // Note that if Kraken paired data (--paired in Kraken2) legacy Illumina read identifier formats
        // with trailing /1 and /2 are stripped of their trails and the reads output does not contain
        // the trails. This mismatches with the input reads (which have /1 and /2) which are therefore
        // not depleted. We do not expect legacy format Illumina reads.

        while let Some(record) = reader.next() {
            let rec = record.map_err(|err| ScrubberError::FastxRecordIO(err))?;
            let rec_id = from_utf8(rec.id()).map_err(|err| ScrubberError::DepletionRecordIdentifier(err))?.split(' ').next().unwrap_or(""); // needletail parses the entire header as identifier (including description)
            
            let to_retain: bool = match extract {
                true => reads.contains(&rec_id.to_string()),
                false => !reads.contains(&rec_id.to_string()) 
            };

            if to_retain {
                rec.write(&mut writer, None).map_err(|err| ScrubberError::FastxRecordIO(err))?;
                read_counts.retained += 1
            } else {
                match extract {
                    true => read_counts.extracted += 1,  // if extraction flag, count extracted, otherwise ...
                    false => read_counts.depleted += 1   // count depleted
                }
            }
            read_counts.total += 1
        }
        Ok(read_counts)
    }
}


// Utility function to get a Needletail reader and Niffler compressed/uncompressed writer
fn get_reader_writer(input: &PathBuf, output: &PathBuf, compression_level: niffler::compression::Level, output_format: Option<niffler::compression::Format>) -> Result<(Box<dyn FastxReader>, Box<dyn std::io::Write>), ScrubberError> {
    // Input output of read files includes compression detection
    let reader = parse_fastx_file(input).map_err(|err| ScrubberError::FastxRecordIO(err))?;

    let file = File::create(&output)?;
    let file_handle = BufWriter::new(file);
    let fmt = match output_format {
        None => niffler::Format::from_path(&output),
        Some(format) => format,
    };

    let writer = niffler::get_writer(Box::new(file_handle), fmt, compression_level).map_err(|err| ScrubberError::DepletionCompressionWriter(err))?;
    
    Ok((reader, writer))
}