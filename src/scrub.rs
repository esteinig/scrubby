use crate::metabuli::MetabuliSeqMode;
use crate::utils::{CompressionExt, FileNameString};
use anyhow::Result;
use chrono::Local;
use clap::crate_version;
use needletail::{parse_fastx_file, FastxReader};
use serde::{Serialize, Serializer};
use std::collections::HashSet;
use std::ffi::OsStr;
use std::fs::File;
use std::fs::{create_dir_all, remove_dir_all};
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
use std::process::Command;
use std::str::from_utf8;
use thiserror::Error;

const JSON_SCHEMA_VERSION: &str = "0.3.0";

#[derive(Error, Debug)]
pub enum ScrubberError {
    /// Represents a failure to execute minimap2
    #[error("failed to detect `minimap2` - is it installed?")]
    MinimapExecutionError,
    /// Represents a failure to run minimap2
    #[error("failed to run `minimap2`")]
    MinimapAlignmentError,
    /// Represents a failure to execute minimap2
    #[error("failed to detect `strobealign` - is it installed?")]
    StrobealignExecutionError,
    /// Represents a failure to run strobealign
    #[error("failed to run `strobealign`")]
    StrobealignAlignmentError,
    /// Represents a failure to execute Kraken2
    #[error("failed to detect `Kraken2` - is it installed?")]
    KrakenExecutionError,
    /// Represents a failure to execute Metabuli
    #[error("failed to detect `Metabuli` - is it installed?")]
    MetabuliExecutionError,
    /// Represents a failure to run Kraken2
    #[error("failed to run `Kraken2`")]
    KrakenClassificationError,
    /// Represents a failure to run Metabuli
    #[error("failed to run `Metabuli`")]
    MetabuliClassificationError,
    /// Represents a failure to count a taxonomic parent during report parsing
    #[error("failed to provide a parent taxon while parsing report from `Kraken2`")]
    KrakenReportTaxonParent,
    /// Represents a failure to convert the read field from string to numeric field in the report file
    #[error("failed to convert the read field in the report from `Kraken2`")]
    KrakenReportReadFieldConversion,
    /// Represents a failure to convert the read field from string to numeric field in the report file
    #[error("failed to convert the direct read field in the report from `Kraken2`")]
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
    /// Represents all other cases of `csv::Error`.
    #[error(transparent)]
    CsvError(#[from] csv::Error),
    /// Indicates that the working directory already exists
    #[error("working directory exists at {0}")]
    WorkdirExists(String),
    /// Indicates a failure to create the working directory
    #[error("working directory path could not be created from {0}")]
    WorkdirCreate(String),
    /// Indicates a failure to obtain an absolute path
    #[error("absolute path could not be obtained from {0}")]
    AbsolutePath(String),
    /// Indicates a failure to parse a file name
    #[error("file name is not valid UTF-8")]
    FileNameInvalidUtf8,
    /// Indicates a failure to parse a file name
    #[error("path does not have a file name")]
    FileNameNotFound,
    /// Indicates failure to parse file with Needletail
    #[error("failed file input/output")]
    FastxRecordIO(#[source] needletail::errors::ParseError),
    /// Indicates failure to obtain compression writer with Niffler
    #[error("failed to get compression writer")]
    DepletionCompressionWriter(#[source] niffler::Error),
    /// Indicates failure to parse record identifier
    #[error("failed to parse record identifier")]
    DepletionRecordIdentifier(#[source] std::str::Utf8Error),
    /// Indicates failure to parse record identifier
    #[error("failed to operate on alignment")]
    ScrubberAlignment(#[source] crate::align::ReadAlignmentError),
    /// Indicates failure with JSON serialization
    #[error("failed to serialize JSON")]
    JsonSerialization(#[source] serde_json::Error),
    /// Indicates a failure to obtain a JSON output file name
    #[error("summary output file name could not be obtained from {0}")]
    JsonNameExtraction(String),
    /// Indicates a failure to pass the correct strobealign mode
    #[error("strobealign mode not supported (0)")]
    StrobealignMode(String),
    /// Indicates a failure to provide the correct strobealign index extension
    #[error("strobealign reference extension must be one of (.fasta|.fa|.sti)")]
    StrobealignReferenceExtension
}

#[derive(Serialize, Clone)]
pub enum ScrubbyTool {
    Kraken2,
    Metabuli,
    Strobealign,
    Minimap2,
    Alignment  // scrub-alignment
}

///
///
///
pub struct Scrubber {
    pub workdir: PathBuf,
    pub json: JsonSummary,
    pub reads: ReadSummary,
    pub output_format: Option<niffler::compression::Format>,
    pub compression_level: niffler::compression::Level,
}

impl Scrubber {
    ///
    ///
    ///
    pub fn new(
        workdir: Option<PathBuf>,
        output_format: Option<niffler::compression::Format>,
        compression_level: niffler::compression::Level,
        settings: Settings,
        force: bool
    ) -> Result<Self, ScrubberError> {
        let _workdir = check_or_create_workdir(workdir, force)?;
        Ok(Self {
            workdir: _workdir,
            json: JsonSummary::new(
                crate_version!().to_string(),
                JSON_SCHEMA_VERSION.to_string(),
                settings,
            ),
            reads: ReadSummary::new(),
            output_format,
            compression_level,
        })
    }
    ///
    ///
    ///
    pub fn test_dependencies(
        &self,
        kraken: bool,
        metabuli: bool,
        minimap2: bool,
        strobealign: bool
    ) -> Result<(), ScrubberError> {

        if kraken {
            // Run the Kraken command
            let output = Command::new("kraken2")
                .args(["-v"])
                .output()
                .map_err(|_| ScrubberError::KrakenExecutionError)?;

            // Ensure command ran successfully
            if !output.status.success() {
                log::error!(
                    "Failed to detect dependency: Kraken2"
                );
                return Err(ScrubberError::KrakenClassificationError);
            }
        }
        if metabuli {
            // Run the Kraken command
            let output = Command::new("metabuli")
                .output()
                .map_err(|_| ScrubberError::KrakenExecutionError)?;

            // Ensure command ran successfully
            if !output.status.success() {
                log::error!(
                    "Failed to detect dependency: Metabuli"
                );
                return Err(ScrubberError::KrakenClassificationError);
            }
        }
        if minimap2 {
            // Run the Kraken command
            let output = Command::new("minimap2")
                .args(["-V"])
                .output()
                .map_err(|_| ScrubberError::KrakenExecutionError)?;

            // Ensure command ran successfully
            if !output.status.success() {
                log::error!(
                    "Failed to detect dependency: minimap2"
                );
                return Err(ScrubberError::KrakenClassificationError);
            }
        }
        if strobealign {
            // Run the Kraken command
            let output = Command::new("strobealign")
                .args(["--version"])
                .output()
                .map_err(|_| ScrubberError::KrakenExecutionError)?;

            // Ensure command ran successfully
            if !output.status.success() {
                log::error!(
                    "Failed to detect dependency: strobealign"
                );
                return Err(ScrubberError::KrakenClassificationError);
            }
        }


        Ok(())
    }
    ///
    ///
    ///
    pub fn run_kraken(
        &self,
        input: &Vec<PathBuf>,
        db_path: &Path,
        db_name: &String,
        db_index: &usize,
        threads: &u32,
        args: &str,
    ) -> Result<(Vec<PathBuf>, String), ScrubberError> {
        // Safely build the arguments for Kraken2
        let kraken_args =
            crate::kraken::get_kraken_command(input, db_path, db_name, db_index, threads, args)?;

        log::info!(
            "Executing taxonomic classification with Kraken2 ({})",
            db_name
        );
        log::debug!("Executing command: {}", &kraken_args.join(" "));

        // Run the Kraken command
        let output = Command::new("kraken2")
            .args(&kraken_args)
            .current_dir(&self.workdir)
            .output()
            .map_err(|_| ScrubberError::KrakenExecutionError)?;

        // Ensure command ran successfully
        if output.status.success() {
            log::info!(
                "Completed taxonomic classification with Kraken2 ({})",
                db_name
            )
        } else {
            log::error!(
                "Failed to run taxonomic classification with Kraken2 ({})",
                db_name
            );
            log::error!("{}", String::from_utf8_lossy(&output.stderr));
            return Err(ScrubberError::KrakenClassificationError);
        }

        let kraken_report = self
            .workdir
            .join(format!("{}-{}.report", db_index, db_name));
        let kraken_reads = self
            .workdir
            .join(format!("{}-{}.kraken", db_index, db_name));

        Ok((Vec::from([kraken_report, kraken_reads]), kraken_args.join(" ")))
    }
    ///
    ///
    ///
    pub fn run_metabuli(
        &self,
        input: &Vec<PathBuf>,
        db_path: &Path,
        db_name: &String,
        db_index: &usize,
        threads: &u32,
        seq_mode: Option<MetabuliSeqMode>,
        args: &str
    ) -> Result<(Vec<PathBuf>, String), ScrubberError> {
        // Safely build the arguments for Metabuli
        let metabuli_args = crate::metabuli::get_metabuli_command(input, db_path, db_name, db_index, threads, seq_mode, args)?;

        log::info!(
            "Executing taxonomic classification with Metabuli ({})",
            db_name
        );
        log::debug!("Executing command: {}", &metabuli_args.join(" "));

        // Run the Kraken command
        let output = Command::new("metabuli")
            .args(&metabuli_args)
            .current_dir(&self.workdir)
            .output()
            .map_err(|_| ScrubberError::MetabuliExecutionError)?;

        // Ensure command ran successfully
        if output.status.success() {
            log::info!(
                "Completed taxonomic classification with Metabuli ({})",
                db_name
            )
        } else {
            log::error!(
                "Failed to run taxonomic classification with Metabuli ({})",
                db_name
            );
            log::error!("{}", String::from_utf8_lossy(&output.stderr));
            return Err(ScrubberError::MetabuliClassificationError);
        }

        let metabuli_report = self
            .workdir
            .join(format!("{}-{}", db_index, db_name))
            .join(format!("{}-{}_report.tsv", db_index, db_name));
        let metabuli_reads = self
            .workdir
            .join(format!("{}-{}", db_index, db_name))
            .join(format!("{}-{}_classifications.tsv", db_index, db_name));

        Ok((Vec::from([metabuli_report, metabuli_reads]), metabuli_args.join(" ")))
    }
    ///
    ///
    ///
    pub fn parse_kraken(
        &self,
        kraken_files: &[PathBuf],
        kraken_taxa: &[String],
        kraken_taxa_direct: &[String],
    ) -> Result<HashSet<String>, ScrubberError> {
        log::info!("Parsing read identifiers from report and classified reads files...");

        let taxids = crate::kraken::get_taxids_from_report(
            kraken_files[0].clone(),
            kraken_taxa,
            kraken_taxa_direct,
        )?;
        crate::kraken::get_taxid_reads_kraken(taxids, kraken_files[1].clone())
    }
    ///
    ///
    ///
    pub fn parse_metabuli(
        &self,
        metabuli_files: &[PathBuf],
        metabuli_taxa: &[String],
        metabuli_taxa_direct: &[String],
    ) -> Result<HashSet<String>, ScrubberError> {
        log::info!("Parsing read identifiers from report and classified reads files...");

        let taxids = crate::kraken::get_taxids_from_report(
            metabuli_files[0].clone(),
            metabuli_taxa,
            metabuli_taxa_direct,
        )?;
        crate::kraken::get_taxid_reads_metabuli(taxids, metabuli_files[1].clone())
    }
    ///
    ///
    ///
    pub fn run_minimap2(
        &self,
        input: &Vec<PathBuf>,
        index_path: &Path,
        index_name: &String,
        index_idx: &usize,
        threads: &u32,
        preset: &String,
        args: &str
    ) -> Result<(PathBuf, String), ScrubberError> {
        let minimap_args = crate::align::get_minimap2_command(
            input, index_path, index_name, index_idx, threads, preset, args
        )?;

        log::info!("Executing read alignment with minimap2 ({})", index_name);
        log::debug!("Executing command: {}", &minimap_args.join(" "));

        let output = Command::new("minimap2")
            .args(&minimap_args)
            .current_dir(&self.workdir)
            .output()
            .map_err(|_| ScrubberError::MinimapExecutionError)?;

        if output.status.success() {
            log::info!("Completed alignment with minimap2 ({})", index_name)
        } else {
            log::error!(
                "Failed to run taxonomic classification with minimap2 ({})",
                index_name
            );
            log::error!("{}", String::from_utf8_lossy(&output.stderr));
            return Err(ScrubberError::MinimapAlignmentError);
        }
        Ok(
            (self
            .workdir
            .join(format!("{}-{}.paf", index_idx, index_name)),
            minimap_args.join(" "))
        )
    }
    ///
    ///
    ///
    pub fn run_strobealign(
        &self,
        input: &Vec<PathBuf>,
        index_path: &Path,
        index_name: &String,
        index_idx: &usize,
        threads: &u32,
        mode: &String,
        args: &str
    ) -> Result<(PathBuf, String), ScrubberError> {
        let strobealign_args = crate::align::get_strobealign_command(
            input, index_path, index_name, index_idx, threads, mode, args
        )?;

        log::info!("Executing read alignment with strobealign ({})", index_name);
        log::debug!("Executing command: {}", &strobealign_args.join(" "));

        let output = Command::new("strobealign")
            .args(&strobealign_args)
            .current_dir(&self.workdir)
            .output()
            .map_err(|_| ScrubberError::StrobealignExecutionError)?;

        // Ensure command ran successfully
        if output.status.success() {
            log::info!("Completed alignment with strobealign ({})", index_name)
        } else {
            log::error!(
                "Failed to run taxonomic classification with strobealign ({})",
                index_name
            );
            log::error!("{}", String::from_utf8_lossy(&output.stderr));
            return Err(ScrubberError::StrobealignAlignmentError);
        }
        match mode.as_str() {
            "map" => Ok((self
                .workdir
                .join(format!("{}-{}.paf", index_idx, index_name)), strobealign_args.join(" "))),
            "align" => Ok((self
                .workdir
                .join(format!("{}-{}.sam", index_idx, index_name)), strobealign_args.join(" "))),
            _ => Err(ScrubberError::StrobealignMode(mode.to_string())),
        }
    }
    ///
    ///
    ///
    pub fn parse_alignment(
        &self,
        alignment: &Path,
        alignment_format: Option<String>,
        min_qaln_len: &u64,
        min_qaln_cov: &f64,
        min_mapq: &u8,
    ) -> Result<HashSet<String>, ScrubberError> {
        log::info!("Parsing read identifiers from alignment file...");
        let alignment = crate::align::ReadAlignment::from(
            alignment,
            *min_qaln_len,
            *min_qaln_cov,
            *min_mapq,
            alignment_format,
        )
        .map_err(ScrubberError::ScrubberAlignment)?;
        Ok(alignment.reads)
    }
    ///
    ///
    ///
    #[allow(clippy::too_many_arguments)]
    pub fn deplete_to_workdir(
        &mut self,
        input: &Vec<PathBuf>,
        reads: &HashSet<String>,
        tool: ScrubbyTool,
        name: &String,
        path: PathBuf,
        idx: &usize,
        extract: &bool,
        command: &str
    ) -> Result<(ReferenceSummary, Vec<PathBuf>), ScrubberError> {
        
        let output = match input.len() {
            2 => Vec::from([
                self.workdir.join(format!("{}-{}_1.fq", idx, name)),
                self.workdir.join(format!("{}-{}_2.fq", idx, name)),
            ]),
            1 => Vec::from([self.workdir.join(format!("{}-{}.fq", idx, name))]),
            _ => return Err(ScrubberError::FileNumberError),
        };

        let mut ref_summary = ReferenceSummary::new(*idx, tool, name.clone(), path, 0, 0, 0, command.to_owned());
        for (i, _) in input.iter().enumerate() {
            let depletor = ReadDepletor::new(self.output_format, self.compression_level)?;
            let read_counts = depletor.deplete(reads, &input[i], &output[i], extract)?;
            ref_summary.files.push(read_counts);
        }
        ref_summary.compute_total();

        Ok((ref_summary, output))
    }
    ///
    ///
    ///
    #[allow(clippy::too_many_arguments)]
    pub fn deplete_to_file(
        &mut self,
        input: &[PathBuf],
        output: &[PathBuf],
        reads: &HashSet<String>,
        tool: ScrubbyTool,
        name: &str,
        path: PathBuf,
        idx: &usize,
        extract: &bool,
        command: &str
    ) -> Result<(ReferenceSummary, Vec<PathBuf>), ScrubberError> {

        let mut ref_summary = ReferenceSummary::new(*idx, tool, name.to_owned(), path, 0, 0, 0, command.to_owned());
        
        for (i, _) in input.iter().enumerate() {
            let depletor = ReadDepletor::new(self.output_format, self.compression_level)?;
            let read_counts = depletor.deplete(reads, &input[i], &output[i], extract)?;
            ref_summary.files.push(read_counts);
        }
        ref_summary.compute_total();

        Ok((ref_summary, output.to_vec()))
    }
    ///
    ///
    ///
    pub fn write_extracted_pipeline_outputs(
        &mut self,
        input_files: Vec<PathBuf>,
        output_files: Vec<PathBuf>,
        reads: &HashSet<String>,
    ) -> Result<(), ScrubberError> {
        let mut total = 0;
        for (i, _) in input_files.iter().enumerate() {
            let depletor = ReadDepletor::new(self.output_format, self.compression_level)?;
            let read_counts = depletor.deplete(reads, &input_files[i], &output_files[i], &true)?;
            total += read_counts.total;
        }
        self.json.update(total);
        Ok(())
    }
    ///
    ///
    ///
    pub fn write_depleted_pipeline_outputs(
        &mut self,
        input_files: Vec<PathBuf>,
        output_files: Vec<PathBuf>,
    ) -> Result<(), ScrubberError> {
        let mut total = 0;
        for (i, _) in input_files.iter().enumerate() {
            // input output are ensured to have same length through command-line interface
           match get_fastx_reader_writer(
                &input_files[i],
                &output_files[i],
                self.compression_level,
                self.output_format,
            )?
            {
                Some(io) => {
                    let (mut reader, mut writer) = io;
                    log::info!(
                        "Writing scrubbed reads to output file: {:?}",
                        output_files[i]
                    );
                    while let Some(record) = reader.next() {
                        let rec = record.map_err(ScrubberError::FastxRecordIO)?;
                        rec.write(&mut writer, None)
                            .map_err(ScrubberError::FastxRecordIO)?;
                        total += 1;
                    }
                },
                None => {
                    log::warn!("All reads seem to have been depleted, writing empty output files");
                }
            };
        }
        self.json.update(total);

        Ok(())
    }
    ///
    ///
    ///
    pub fn write_summary(&self, json: Option<PathBuf>) -> Result<(), ScrubberError> {
        if let Some(json_file) = json {
            match json_file.file_name() {
                Some(name) => match name == OsStr::new("-") {
                    true => self.print_json()?,
                    false => self.write_json(json_file)?,
                },
                None => {
                    return Err(ScrubberError::JsonNameExtraction(format!(
                        "{:?}",
                        json_file
                    )))
                }
            };
        }

        Ok(())
    }
    ///
    /// 
    /// 
    pub fn write_read_summary(&self, tsv: Option<PathBuf>) -> Result<(), ScrubberError> {

        if let Some(tsv_file) = tsv {
            log::info!("Writing read summary to: {}", tsv_file.display());
            let mut wrtr = csv::WriterBuilder::new().delimiter(b'\t').has_headers(true).from_path(tsv_file)?;
            for record in &self.reads.records {
                wrtr.serialize(record)?;
            }
        }
        
        Ok(())

    }
    ///
    ///
    ///
    pub fn write_json(&self, output: PathBuf) -> Result<(), ScrubberError> {
        log::info!("Writing JSON summary to: {}", output.display());
        let mut file = File::create(&output)?;
        let json_string =
            serde_json::to_string_pretty(&self.json).map_err(ScrubberError::JsonSerialization)?;
        write!(file, "{}", json_string)?;
        Ok(())
    }
    ///
    ///
    ///
    pub fn print_json(&self) -> Result<(), ScrubberError> {
        let json_string =
            serde_json::to_string_pretty(&self.json).map_err(ScrubberError::JsonSerialization)?;
        println!("{}", json_string);
        Ok(())
    }
    ///
    ///
    ///
    pub fn clean_up(&self, keep: bool) -> Result<(), ScrubberError> {
        match keep {
            true => {
                log::info!("Keeping working directory and intermediary files");
            }
            false => {
                log::info!("Deleting working directory and intermediary files");
                remove_dir_all(&self.workdir)?;
            }
        }
        Ok(())
    }
}

/// Checks if work directory exists and otherwise creates it
pub fn check_or_create_workdir(workdir: Option<PathBuf>, force: bool) -> Result<PathBuf, ScrubberError> {
    let _workdir = match workdir {
        Some(path) => path,
        None => PathBuf::from(format!("Scrubby-{}", Local::now().format("%Y%m%dT%H%M%S"))),
    };
    let workdir_msg = format!("{:?}", _workdir);
    if !&_workdir.exists() {
        create_dir_all(&_workdir).map_err(|_| ScrubberError::WorkdirCreate(workdir_msg))?;
        let abs_workdir = _workdir
            .canonicalize()
            .map_err(|_| ScrubberError::AbsolutePath(format!("{:?}", _workdir)))?;
        Ok(abs_workdir)
    } else {
        if force {
            log::info!("Force overwrite active, removing existing directory");
            remove_dir_all(&_workdir)?;
            create_dir_all(&_workdir).map_err(|_| ScrubberError::WorkdirCreate(workdir_msg))?;
            let abs_workdir = _workdir
                .canonicalize()
                .map_err(|_| ScrubberError::AbsolutePath(format!("{:?}", _workdir)))?;
            Ok(abs_workdir)
        } else {
            Err(ScrubberError::WorkdirExists(workdir_msg))
        }
    }
}

/*
=========================
Read depletion/extraction
=========================
*/

// Struct to hold the read depletion for
// each reference/database
#[derive(Serialize)]
pub struct JsonSummary {
    pub version: String,
    pub schema_version: String,
    pub summary: Summary,
    pub settings: Settings,
    pub pipeline: Vec<ReferenceSummary>,
}

impl JsonSummary {
    ///
    ///
    ///
    pub fn new(version: String, schema_version: String, settings: Settings) -> Self {
        Self {
            version,
            schema_version,
            summary: Summary {
                total: 0,
                depleted: 0,
                extracted: 0,
            },
            settings,
            pipeline: Vec::new(),
        }
    }
    ///
    ///
    ///
    pub fn update(&mut self, total: u64) {
        // Total reads are passed during the output iteration
        // because we don't always have them available from the
        // pipeline (e.g. in the individual scrubbing tasks)
        match self.settings.extract {
            true => {
                self.summary.total = total;
            }
            false => match self.pipeline.get(0) {
                Some(first_summary) => self.summary.total = first_summary.total,
                None => self.summary.total = total,
            },
        }

        for summary in self.pipeline.iter() {
            self.summary.depleted += summary.depleted;
            self.summary.extracted += summary.extracted;
        }
    }
}

// Struct to hold a compiled run summary
#[derive(Serialize)]
pub struct Summary {
    pub total: u64,
    pub depleted: u64,
    pub extracted: u64,
}

// Struct to hold the provided settings
#[derive(Serialize)]
pub struct Settings {
    pub taxa: Vec<String>,
    pub taxa_direct: Vec<String>,
    pub min_len: u64,
    pub min_cov: f64,
    pub min_mapq: u8,
    pub extract: bool,
}
impl Settings {
    pub fn new(
        taxa: Vec<String>,
        taxa_direct: Vec<String>,
        min_len: u64,
        min_cov: f64,
        min_mapq: u8,
        extract: bool,
    ) -> Self {
        Self {
            taxa,
            taxa_direct,
            min_len,
            min_cov,
            min_mapq,
            extract,
        }
    }
}

// A summary struct to hold counts
// for the read depletion/extraction
#[derive(Clone, Serialize)]
pub struct FileSummary {
    pub total: u64,
    pub depleted: u64,
    pub extracted: u64,
    #[serde(serialize_with = "serialize_path_as_filename")]
    pub input_file: PathBuf,
    #[serde(serialize_with = "serialize_path_as_filename")]
    pub output_file: PathBuf,
}

// Struct to hold the read depletion for
// each reference/database
#[derive(Clone, Serialize)]
pub struct ReferenceSummary {
    pub index: usize,
    pub tool: ScrubbyTool,
    pub name: String,
    #[serde(serialize_with = "serialize_path_as_filename")]
    pub path: PathBuf,
    pub total: u64,
    pub depleted: u64,
    pub extracted: u64,
    pub files: Vec<FileSummary>,
    pub command: String
}

impl ReferenceSummary {
    ///
    ///
    ///
    pub fn new(
        index: usize,
        tool: ScrubbyTool,
        name: String,
        path: PathBuf,
        total: u64,
        depleted: u64,
        extracted: u64,
        command: String
    ) -> Self {
        Self {
            index,
            tool,
            name,
            path,
            total,
            depleted,
            extracted,
            files: Vec::new(),
            command: command
        }
    }
    ///
    ///
    ///
    pub fn compute_total(&mut self) {
        for file_summary in self.files.iter() {
            self.total += file_summary.total;
            self.depleted += file_summary.depleted;
            self.extracted += file_summary.extracted;
        }
    }
}
// Custom serializer for outputting only file names not full paths
fn serialize_path_as_filename<S>(path: &PathBuf, serializer: S) -> Result<S::Ok, S::Error>
where
    S: Serializer,
{
    let filename = PathBuf::file_name_string(path).map_err(serde::ser::Error::custom)?;
    filename.serialize(serializer)
}


// Struct to hold the read depletion for
// each reference/database
#[derive(Clone, Serialize)]
pub struct ReadSummary {
    pub records: Vec<ReadSummaryRecord>
}
impl ReadSummary {
    pub fn new() -> Self {
        Self {
            records: Vec::new()
        }
    }
    pub fn add(&mut self, reads: &HashSet<String>, tool: ScrubbyTool, db: &str) -> () {
        for read in reads {
            self.records.push(ReadSummaryRecord { id: read.clone(), tool: tool.clone(), db: db.to_string() })
        }
    }
}

// Struct to hold the read depletion for
// each reference/database
#[derive(Clone, Serialize)]
pub struct ReadSummaryRecord {
    pub id: String,
    pub tool: ScrubbyTool,
    pub db: String
}

/// Read depletion struct
///
/// This struct can be initialised with `niffler ::compression::Format` which
/// optionally specifies the output format and a `niffler::compression::Level`
/// which specified the output compression level.
pub struct ReadDepletor {
    output_format: Option<niffler::compression::Format>,
    compression_level: niffler::compression::Level,
}

impl ReadDepletor {
    ///
    ///
    ///
    pub fn new(
        output_format: Option<niffler::compression::Format>,
        compression_level: niffler::compression::Level,
    ) -> Result<Self, ScrubberError> {
        Ok(Self {
            output_format,
            compression_level,
        })
    }
    /// Deplete reads from an input read file
    ///
    /// This method depletes (`extract = false`) or extracts (`extract = true`) an input
    /// read file which may be compressed. It checks if the read identifier is contained
    /// within the `HashSet`of read identifiers provided. It also counts the total,
    /// depleted or extracted, and the retained reads and returns a `ReadCounts` object
    /// which can be added to the `ReadCountSummary` to be output to JSON.
    ///
    pub fn deplete(
        &self,
        reads: &HashSet<String>,
        input: &PathBuf,
        output: &PathBuf,
        extract: &bool,
    ) -> Result<FileSummary, ScrubberError> {

        let mut read_counts = FileSummary {
            total: 0,
            depleted: 0,
            extracted: 0,
            input_file: input.to_path_buf(),
            output_file: output.to_path_buf(),
        };

        // Input output of read files includes compression detection
        let (mut reader, mut writer) = match get_fastx_reader_writer(input, output, self.compression_level, self.output_format)? {
            Some(io) => io,
            None => {
                log::warn!("Failed to deplete input file, returning zero counts");
                return Ok(read_counts)
            }
        };

        while let Some(record) = reader.next() {
            let rec = record.map_err(ScrubberError::FastxRecordIO)?;
            
            let rec_id = from_utf8(rec.id())
                .map_err(ScrubberError::DepletionRecordIdentifier)?
                .split(' ')
                .next()
                .unwrap_or(""); // needletail parses the entire header as identifier (including description)

            let to_retain: bool = match extract {
                true => reads.contains(&rec_id.to_string()),
                false => !reads.contains(&rec_id.to_string()),
            };

            if to_retain {
                rec.write(&mut writer, None)
                    .map_err(ScrubberError::FastxRecordIO)?;
                if *extract {
                    read_counts.extracted += 1
                }
            } else {
                read_counts.depleted += 1
            }
            read_counts.total += 1
        }
        Ok(read_counts)
    }
}

#[allow(clippy::type_complexity)]
// Utility function to get a Needletail reader and Niffler compressed/uncompressed writer
fn get_fastx_reader_writer(
    input: &PathBuf,
    output: &PathBuf,
    compression_level: niffler::compression::Level,
    output_format: Option<niffler::compression::Format>,
) -> Result<Option<(Box<dyn FastxReader>, Box<dyn std::io::Write>)>, ScrubberError> {

    // Input output of read files includes compression detection
    let reader = parse_fastx_file(input).ok();

    match reader {
        Some(reader) => {

            let file = File::create(output)?;
            let file_handle = BufWriter::new(file);
            let fmt = match output_format {
                None => niffler::Format::from_path(output),
                Some(format) => format,
            };
        
            let writer = niffler::get_writer(Box::new(file_handle), fmt, compression_level)
                .map_err(ScrubberError::DepletionCompressionWriter)?;
        
            Ok(Some((reader, writer)))
        },
        None => {
            log::info!("Could not parse input file, it may be empty: {}", &input.display());
            let _ = File::create(output)?; // Creates the empty output file regardless! Important for continuity of pipeline
            Ok(None)
        }
    }

}
