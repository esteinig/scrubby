use anyhow::{Result, Context};
use std::{path::PathBuf, process::Output, collections::HashSet, io::{BufRead, BufReader}, fs::File};
use needletail::{parse_fastx_file, Sequence};
use std::str::from_utf8;
use std::process::Command;
use thiserror::Error;
use chrono::Local;
use std::fs::{create_dir_all};

// See: https://nick.groenen.me/posts/rust-error-handling/

/// WordCountError enumerates all possible errors returned by this library.
#[derive(Error, Debug)]
pub enum ScrubberError {
    /// Represents a failure to execute Kraken2
    #[error("failed to run Kraken2 - is it installed?")]
    KrakenExecutionError,
    /// Represents a failure to run Kraken2
    #[error("failed to run Kraken2")]
    KrakenClassificationError,
    /// Represents a failure in the correct length of the input file vector
    #[error("input file error - incorrect number of input files")]
    FileNumberError,
    /// Represents a failure to convert a PathBuf converted OsStr into a String
    /// This is because into_string() returns Result<String, OsString)
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
    /// Indicates a failure to obtain an absolute path
    #[error("database name could not be obtained from {0}")]
    DatabaseName(String)
}

pub struct Scrubber {
    workdir: PathBuf
}

impl Scrubber {
    pub fn new(workdir: Option<PathBuf>) -> Result<Self, ScrubberError> {
        let _workdir = check_or_create_workdir(workdir)?;
        Ok(Self { workdir: _workdir })
    }
    pub fn run_kraken(
        &self,
        input: &Vec<PathBuf>,
        db: PathBuf,
        threads: u32,
    ) -> Result<(), ScrubberError>{
        
        // Kraken2 database name
        let db_name = match db.file_name(){
            Some(name) => name.to_os_string().into_string().map_err(|_| ScrubberError::InvalidFilePathConversion)?,
            None => return Err(ScrubberError::DatabaseName(format!("{:?}", db)))
        };

        // Safely build the arguments for Kraken2
        let kraken_args = get_kraken_command(input, db, threads, &db_name)?;

        // Run the Kraken command
        let output = Command::new("kraken2")
            .args(kraken_args)
            .current_dir(&self.workdir)
            .output()
            .map_err(|_| ScrubberError::KrakenExecutionError)?;

        // Ensure command ran successfully
        if output.status.success() {
            log::info!("Kraken2 - completed taxonomic assignment ({})", db_name)
        } else {
            log::info!("Kraken2 - failed to run taxonomic assignment ({})", db_name);
            let err_msg = get_kraken_err_msg(output)?;
            log::info!("Error from {}", err_msg);
            return Err(ScrubberError::KrakenClassificationError)
        }
        
        Ok(())
    }
    pub fn deplete_kraken(
        &self,
        input: &Vec<PathBuf>,
        kraken_taxa: &Vec<String>,
        kraken_taxa_direct: &Vec<String>
    ) -> Result<(), ScrubberError>{

        log::info!("Kraken2 - depleting reads across selected taxa.");

        Ok(())
    }
}



/// Checks if work directory exists and otherwise creates it
///  
/// # Errors
/// A [`ScrubberError::WorkdirExists`](#clierror) is returned if the directory already exists
/// A [`ScrubberError::WorkdirCreate`](#clierror) is returned if the directory cannot be created
/// A [`ScrubberError::AbsolutePath`](#clierror) is returned if the directory path cannot be canonicalized
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

/// Builds the Kraken2 command from the input configuration
/// 
/// # Errors
/// A [`ScrubberError::InvalidFilePathConversion`](#clierror) is returned if one of the input paths could not be converted to a string
/// A [`ScrubberError::FileNumberError`](#clierror) is returned if for some arcane reason the input file vector is not the correct length
pub fn get_kraken_command(input: &Vec<PathBuf>, db: PathBuf, threads: u32, db_name: &str) -> Result<Vec<String>, ScrubberError> {
    
    let kraken_db_path = db.into_os_string().into_string().map_err(|_| ScrubberError::InvalidFilePathConversion)?;
    let kraken_threads_arg = threads.to_string();

    let file_arg = match input.len() {
        2 => {
            let file1 = input[0].clone().into_os_string().into_string().map_err(|_| ScrubberError::InvalidFilePathConversion)?;
            let file2 = input[1].clone().into_os_string().into_string().map_err(|_| ScrubberError::InvalidFilePathConversion)?;
            [Some(file1), Some(file2)]
        },
        1 => [Some(input[0].clone().into_os_string().into_string().map_err(|_| ScrubberError::InvalidFilePathConversion)?), None],
        _ => return Err(ScrubberError::FileNumberError),
    };
    let paired_arg = match input.len() {
        2 => Some("--paired"),
        1 => None,
        _ => return Err(ScrubberError::FileNumberError),
    };
   
    let mut kraken_args = Vec::from([
        "--threads".to_string(),
        kraken_threads_arg, 
        "--db".to_string(), 
        kraken_db_path,
        "--output".to_string(),
        format!("{}.kraken", db_name),
        "--report".to_string(),
        format!("{}.report", db_name)
    ]);

    match paired_arg {
        Some(value) => kraken_args.push(value.to_string()),
        None => {}
    };
    for file in file_arg {
        match file {
            Some(value) => kraken_args.push(value),
            None => {}
        }
    };

    Ok(kraken_args)
}

/// Parses the error message from Kraken2
///  
/// # Errors
/// A [`ScrubberError::KrakenClassificationError`](#clierror) is returned if parsing the error message fails
pub fn get_kraken_err_msg(cmd_output: Output) -> Result<String, ScrubberError>{
    let err_out = String::from_utf8_lossy(&cmd_output.stderr);
    match err_out.lines().nth(0){
        Some(msg) => Ok(msg.to_string()),
        None => return Err(ScrubberError::KrakenClassificationError)
    }
}

enum TaxonomicLevel {
    Unclassified,
    Root,
    Domain,
    Phylum,
    Class,
    Order,
    Family,
    Genus,
    Species
}


pub fn parse_kraken_taxids(
    // Kraken taxonomic report
    kraken_report: PathBuf,
    // Kraken read classifications
    kraken_reads: PathBuf,
    kraken_taxa: Vec<String>,
    kraken_taxa_direct: Vec<String>
) -> Result<HashSet<String>, ScrubberError>{

    let report = BufReader::new(File::open(kraken_report)?);

    // Make sure no trailign whitespaces are input by user 
    let kraken_taxa: Vec<String> = kraken_taxa.into_iter().map(|x| x.trim().to_string()).collect();
    let kraken_taxa_direct: Vec<String> = kraken_taxa_direct.into_iter().map(|x| x.trim().to_string()).collect();

    let mut taxids: HashSet<String> = HashSet::new();
    let mut extract_taxids: bool = false;

    'report: for line in report.lines() {  
        // Iterate over the lines in the report file - it is not known which domain comes
        // first, so we set a flag when we want to extract identifiers, until the next taxonomic
        // level report line of the same kind or above, when the requested domains are checked again
        let record: KrakenReportRecord = KrakenReportRecord::from_str(line?)?;
        
        // Skip all records that do not have reads directly assigned to it!
        if record.reads_direct == 0 {
            continue 'report;
        }

        // Add the direct taxon identifier if the record matches by taxonomic name or identifer 
        // (and has reads assigned directly)
        if kraken_taxa_direct.contains(&record.tax_name) || kraken_taxa_direct.contains(&record.tax_id) {
            taxids.insert(record.tax_id);
        }

        match record.tax_level.as_str() {

            "D" => {
                if domains_to_deplete.contains(&record.tax_name) {
                    extract_taxids = true;
                    taxids.insert(record.tax_id);
                } else {
                    extract_taxids = false;
                }
            },
            &_ => {
                if extract_taxids {
                    taxids.insert(record.tax_id);
                }
            }
        }
    }


    // Extraction of read identifiers extracted from the report or added directly above
    let file = BufReader::new(File::open(&path)?);
    let mut reads: HashSet<String> = HashSet::new();
    for line in file.lines(){            
        let record: KrakenReadRecord = KrakenReadRecord::from_str(line?)?;
        if taxids.contains(&record.tax_id){
            reads.insert(record.read_id.clone());
        }
    }

    println!("Depleting a total of {:} reads", reads.len());

    Ok(reads)

}

/*
==============
Record Structs
==============
*/

/// Kraken read classification record - we are handling taxonomic identifiers
/// as strings in case they are not numeric (e.g. GTDB)
#[derive(Debug, Clone)]
pub struct KrakenReadRecord {
    pub classified: bool,
    pub read_id: String,
    pub tax_id: String,
    pub read_len: String,
    pub annotation: String,
}

impl KrakenReadRecord {
    pub fn from_str(kraken_line: String) -> Result<Self, ScrubberError> {
        let fields: Vec<&str> = kraken_line.split('\t').collect();
        
        let _classified = match fields[0] {
            "U" => false,
            "C" => true,
            _ => false
        };

        let record = Self {
            classified: _classified,
            read_id: fields[1].trim().to_string(),
            tax_id: fields[2].trim().to_string(),
            read_len: fields[3].trim().to_string(),
            annotation: fields[4].trim().to_string()
        };

        Ok(record)
    }
}

/// Kraken read classification record
#[derive(Debug, Clone)]
pub struct KrakenReportRecord {
    pub fraction: String,
    pub reads: u64,
    pub reads_direct: u64,
    pub tax_level: String,
    pub tax_id: String,
    pub tax_name: String,
}

impl KrakenReportRecord {
    // Create a record from a parsed line
    pub fn from_str(report_line: String) -> Result<Self, ScrubberError> {
        let fields: Vec<&str> = report_line.split('\t').collect();

        let record = Self {
            fraction: fields[0].to_string(),
            reads: fields[1].parse::<u64>()?,
            reads_direct: fields[2].parse::<u64>()?,
            tax_level: fields[3].trim().to_string(),
            tax_id: fields[4].trim().to_string(),
            tax_name: fields[5].trim().to_string()
        };

        Ok(record)
    }
}