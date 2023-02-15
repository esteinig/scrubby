
use std::fmt;
use chrono::Local;
use std::fs::File;
use anyhow::Result;
use serde::Serialize;
use thiserror::Error;
use std::path::PathBuf;
use std::str::from_utf8;
use std::process::Output;
use std::process::Command;
use std::fs::create_dir_all;
use std::collections::HashSet; 
use needletail::parse_fastx_file;
use crate::utils::CompressionExt;
use std::io::{BufRead, BufReader, BufWriter};


#[derive(Error, Debug)]
pub enum ScrubberError {
    /// Represents a failure to execute Kraken2
    #[error("failed to run Kraken2 - is it installed?")]
    KrakenExecutionError,
    /// Represents a failure to run Kraken2
    #[error("failed to run Kraken2")]
    KrakenClassificationError,
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
    DepletionFastxParser(#[source] needletail::errors::ParseError),
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
        let kraken_args = get_kraken_command(input, db_path, db_name, db_index, threads)?;

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
            let err_msg = get_kraken_err_msg(output)?;
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

        let msg_word = match extract { true => "Extracting", false => "Depleting" };
        log::info!("{} classified reads from Kraken2", &msg_word);
        
        let reads = parse_kraken_files(
            kraken_files[0].clone(),
            kraken_files[1].clone(),
            kraken_taxa,
            kraken_taxa_direct
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
        for (i, input_file) in input.iter().enumerate() {

            log::info!("{} reads {:#?} into {:#?}", &msg_word, &input_file, &output[i]);

            // Initiate the depletion operator and deplete/extract the reads identifiers parsed from Kraken
            let depletor = ReadDepletion::new(self.output_format, self.compression_level)?;
            let read_counts = depletor.deplete(&reads, &input[i], &output[i], extract)?;

            read_summary.reads.push(read_counts);
        };

        Ok(output)
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

/// Builds the Kraken2 command from the input configuration
/// 
/// # Errors
/// A [`ScrubberError::InvalidFilePathConversion`](#scrubbererror) is returned if one of the input paths could not be converted to a string
/// A [`ScrubberError::FileNumberError`](#scrubbererror) is returned if for some arcane reason the input file vector is not the correct length
pub fn get_kraken_command(input: &Vec<PathBuf>, db_path: &PathBuf, db_name: &str, db_index: &usize, threads: &u32) -> Result<Vec<String>, ScrubberError> {
    
    let kraken_db_path = db_path.to_path_buf().into_os_string().into_string().map_err(|_| ScrubberError::InvalidFilePathConversion)?;
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
        format!("{}-{}.kraken", db_index, db_name),
        "--report".to_string(),
        format!("{}-{}.report", db_index, db_name)
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
/// A [`ScrubberError::KrakenClassificationError`](#scrubbererror) is returned if parsing the error message fails
pub fn get_kraken_err_msg(cmd_output: Output) -> Result<String, ScrubberError>{
    let err_out = String::from_utf8_lossy(&cmd_output.stderr);
    match err_out.lines().nth(0){
        Some(msg) => Ok(msg.to_string()),
        None => return Err(ScrubberError::KrakenClassificationError)
    }
}
/// Taxonomic level enumeration
///
/// Provides an integer value for comparison of levels
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
pub enum TaxonomicLevel {
    None,
    Unspecified,
    Unclassified,
    Root,
    Domain,
    Kingdom,
    Phylum,
    Class,
    Order,
    Family,
    Genus,
    Species
}
impl fmt::Display for TaxonomicLevel {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            TaxonomicLevel::None => write!(f, "None"),
            TaxonomicLevel::Unspecified => write!(f, "Unspecified"),
            TaxonomicLevel::Unclassified => write!(f, "Unclassified"),
            TaxonomicLevel::Root => write!(f, "Root"),
            TaxonomicLevel::Domain => write!(f, "Domain"),
            TaxonomicLevel::Kingdom => write!(f, "Kingdom"),
            TaxonomicLevel::Phylum => write!(f, "Phylum"),
            TaxonomicLevel::Class => write!(f, "Class"),
            TaxonomicLevel::Order => write!(f, "Order"),
            TaxonomicLevel::Family => write!(f, "Family"),
            TaxonomicLevel::Genus => write!(f, "Genus"),
            TaxonomicLevel::Species => write!(f, "Species"),
        }
    }
}

/// Parse the Kraken output report and read file
/// 
/// This functions implements the logic to first parse the record file, and extract any
/// directly specified taxon by name or identifier. It also allows for extraction of any
/// taxon sub-leves of taxa by name or identifier, for example when specifying the
/// domain `Eukaryota` it will parse all sub-levels of `Eukaryota` until the next full 
/// domain level is reached (which prevents sub-specifications of domains to be excluded,
/// such as D1, D2, D3)
/// 
/// It should be ensured that the underlying database taxonomy does not incorrectly specify
/// sub-levels as the triggering level - for example, when specifying `Eukaryota` the 
/// 16S rRNA SILVA database incorrectly specifies `Holozoa` at the domain level, so it
/// should be includedin the taxa to dpeletye to ensure all sequences below `Eukaryota`'
/// and `Holozoa` are excluded.
/// 
/// Only taxa with directly assigned reads are included in the to-deplete list.
/// 
/// Returns a `HashSet` of read identifiers parsed from the read classification file
/// where the classification must be included in the to-deplete list derived from the report file.
/// 
/// # Errors
/// A [`ScrubberError::KrakenReportReadFieldConversion`](#scrubbererror) is returned if the read field in the report file cannot be converted into `u64`
/// A [`ScrubberError::KrakenReportDirectReadFieldConversion`](#scrubbererror) is returned if the direct read field in the report file cannot be converted into `u64`
pub fn parse_kraken_files(
    // Kraken taxonomic report
    kraken_report: PathBuf,
    // Kraken read classifications
    kraken_reads: PathBuf,
    kraken_taxa: &Vec<String>,
    kraken_taxa_direct: &Vec<String>
) -> Result<HashSet<String>, ScrubberError>{

    let report = BufReader::new(File::open(kraken_report)?);

    // Make sure no trailign whitespaces are input by user - these are taxon names or taxon identifiers to deplete
    let kraken_taxa: Vec<String> = kraken_taxa.into_iter().map(|x| x.trim().to_string()).collect();
    let kraken_taxa_direct: Vec<String> = kraken_taxa_direct.into_iter().map(|x| x.trim().to_string()).collect();

    let mut taxids: HashSet<String> = HashSet::new();
    let mut extract_taxlevel: TaxonomicLevel = TaxonomicLevel::None; // make sure this makes sense to initialize 

    'report: for line in report.lines() {  
        // Iterate over the lines in the report file - it is not known which domain comes
        // first, so we set a flag when we want to extract identifiers, until the next taxonomic
        // level report line of the same kind or above, when the requested domains are checked again
        let record: KrakenReportRecord = KrakenReportRecord::from_str(line?)?;

        let tax_level = get_tax_level(&record);

        // Add the direct taxon identifier if the record matches by taxonomic name or identifer 
        // (and has reads assigned directly) - this is always the case when the direct taxon is found
        if kraken_taxa_direct.contains(&record.tax_name) || kraken_taxa_direct.contains(&record.tax_id) {
            log::info!("Detected direct taxon to deplete ({} : {} : {} : {} : {})", &tax_level.to_string(), &record.tax_level, &record.tax_id, &record.tax_name, &record.reads_direct);
            taxids.insert(record.tax_id.clone());
        }

        // If taxon level is above Domain - do not allow it to be processed with the block statement below
        // this is to prevent failure of the logic implemented to parse the report sub-levels of a given
        // taxonomic name or identifier

        if tax_level < TaxonomicLevel::Domain { // Unspecified, Unclassified, Root --> all should be given directly!
            log::warn!("Detected taxon level below `Domain` - ignored in sublevel depletion ({} : {} : {} : {} : {})", &tax_level.to_string(), &record.tax_level, &record.tax_id, &record.tax_name, &record.reads_direct);
            continue 'report;
        } 

        if kraken_taxa.contains(&record.tax_name) || kraken_taxa.contains(&record.tax_id) {
            log::info!("Detected taxon level to deplete ({} : {} : {} : {} : {})", &tax_level.to_string(), &record.tax_level, &record.tax_id, &record.tax_name, &record.reads_direct);
            // If the current record is in the vector of taxa (and following sub-taxa) to deplete, switch on the extraction flag 
            // and set the current tax level as a flag for stopping the extraction in subsequent records that are below or equal
            // to this tax level
            extract_taxlevel = tax_level;
            log::debug!("Setting taxon level for extraction of sub-levels to {} ({})", extract_taxlevel.to_string(), &record.tax_name);
            // Skip all records that do not have reads directly assigned to it!
            if record.reads_direct > 0 {
                taxids.insert(record.tax_id);
            }
        } else {
            if extract_taxlevel == TaxonomicLevel::None { // guard against no taxa given on initial loop
                log::debug!("Ignoring record ({} : {} : {} : {} : {})", &tax_level.to_string(), &record.tax_level, &record.tax_id, &record.tax_name, &record.reads_direct);
                continue 'report;
            }
            // If the current record is not the depletion list, first check if the taxon level indicates we need to stop - this
            // is the case if the tax level is the same or below the extraction flagged tax level set when a record was found to
            // start depletion 
            if (tax_level <= extract_taxlevel) && (record.tax_level.len() == 1) { //  guard against premature sub-level reset (e.g. D1, D2, D3)
                // Unset the extraction flag and reset the taxonomic level
                // to the lowest setting (None - does not ocurr in report)
                log::debug!("Detected taxon level for sub-level reset ({} : {} : {} : {} : {})", &tax_level.to_string(), &record.tax_level, &record.tax_id, &record.tax_name, &record.reads_direct);
                extract_taxlevel = TaxonomicLevel::None;
            } else {
                // Otherwise the taxonomic level is below the one set in the flag and the taxon should be depleted
                // Skip all records that do not have reads directly assigned to it!
                if record.reads_direct > 0 {
                    log::debug!("Detected taxon sub-level with reads to deplete ({} : {} : {} : {} : {})", &tax_level.to_string(), &record.tax_level, &record.tax_id, &record.tax_name, &record.reads_direct);
                    taxids.insert(record.tax_id);
                }
            }

        }
    }
    // HashSet of read identifiers for later depletion
    let mut reads: HashSet<String> = HashSet::new();

    // Extraction of read identifiers extracted from the report or added directly above
    let file = BufReader::new(File::open(&kraken_reads)?);
    for line in file.lines(){            
        let record: KrakenReadRecord = KrakenReadRecord::from_str(line?)?;
        if taxids.contains(&record.tax_id){
            reads.insert(record.read_id.clone());
        }
    }

    log::info!("Parsed Kraken output files; a total of {} reads will be depleted", reads.len());

    Ok(reads)

}

/// A utility function to extract a non-specific tax level into an enumeration
fn get_tax_level(record: &KrakenReportRecord) -> TaxonomicLevel {
    let tax_level_str = &record.tax_level;

    if tax_level_str.starts_with("U")       { return TaxonomicLevel::Unclassified }
    else if tax_level_str.starts_with("R")  { return TaxonomicLevel::Root }
    else if tax_level_str.starts_with("D")  { return TaxonomicLevel::Domain }
    else if tax_level_str.starts_with("K")  { return TaxonomicLevel::Kingdom }
    else if tax_level_str.starts_with("P")  { return TaxonomicLevel::Phylum }
    else if tax_level_str.starts_with("C")  { return TaxonomicLevel::Class }
    else if tax_level_str.starts_with("O")  { return TaxonomicLevel::Order }
    else if tax_level_str.starts_with("F")  { return TaxonomicLevel::Family }
    else if tax_level_str.starts_with("G")  { return TaxonomicLevel::Genus }
    else if tax_level_str.starts_with("S")  { return TaxonomicLevel::Species }
    else                                    { return TaxonomicLevel::Unspecified }
}

/*
==============
Record Structs
==============
*/

// Kraken read classification record - we are handling taxonomic identifiers
// as strings in case they are not numeric (e.g. GTDB)
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

// Kraken read classification record
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
            reads: fields[1].parse::<u64>().map_err(|_| ScrubberError::KrakenReportReadFieldConversion)?,
            reads_direct: fields[2].parse::<u64>().map_err(|_| ScrubberError::KrakenReportDirectReadFieldConversion)?,
            tax_level: fields[3].trim().to_string(),
            tax_id: fields[4].trim().to_string(),
            tax_name: fields[5].trim().to_string()
        };

        Ok(record)
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
        let mut reader = parse_fastx_file(input).map_err(|err| ScrubberError::DepletionFastxParser(err))?;

        let file = File::create(&output)?;
        let file_handle = BufWriter::new(file);
        let fmt = match self.output_format {
            None => niffler::Format::from_path(&output),
            Some(format) => format,
        };

        let mut writer_retained = niffler::get_writer(Box::new(file_handle), fmt, self.compression_level).map_err(|err| ScrubberError::DepletionCompressionWriter(err))?;
        

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
            let rec = record.map_err(|err| ScrubberError::DepletionFastxParser(err))?;
            let rec_id = from_utf8(rec.id()).map_err(|err| ScrubberError::DepletionRecordIdentifier(err))?.split(' ').next().unwrap_or(""); // needletail parses the entire header as identifier (including description)
            
            let to_retain: bool = match extract {
                true => reads.contains(&rec_id.to_string()),
                false => !reads.contains(&rec_id.to_string()) 
            };

            if to_retain {
                rec.write(&mut writer_retained, None).map_err(|err| ScrubberError::DepletionFastxParser(err))?;
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
