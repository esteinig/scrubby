use std::path::PathBuf;
use clap::{crate_version, ArgGroup, Args, Parser, Subcommand};

use crate::prelude::*;
use crate::error::ScrubbyError;

#[derive(Debug, Parser)]
#[command(
    author="Eike Steinig (@esteinig)", 
    version=crate_version!(), 
    about="Taxonomic read depletion for clinical metagenomic diagnostics",
    help_template="\
{before-help}{name} {version} 
{author-with-newline}
{about-with-newline}
{usage-heading} {usage}

{all-args}{after-help}
"
)]
#[command(styles=get_styles())]
#[command(arg_required_else_help(true))]
pub struct App {
    /// Output logs to file instead of terminal
    ///
    /// Specify a file path to store the output logs. If not provided,
    /// logs will be displayed in the terminal. 
    #[arg(short, long)]
    pub log_file: Option<PathBuf>,

    #[clap(subcommand)]
    pub command: Commands,
}

/// Enumeration of available commands for Scrubby.
#[derive(Debug, Subcommand)]
pub enum Commands {
    /// Deplete or extract reads using aligners or k-mer classifiers.
    Clean(CleanArgs),
    /// Deplete or extract reads from k-mer classifier outputs (Kraken2, Metabuli).
    Classifer(ClassifierArgs),
    /// Deplete or extract reads from aligner output with additional filters (SAM/BAM/PAF).
    Alignment(AlignmentArgs),
    /// List available indices and download files for aligners and classfiers.
    Download(DownloadArgs),
}

/// Command-line arguments for the cleaning operation
///
/// This struct defines the arguments required for performing the cleaning 
/// operation using Scrubby. It includes options for specifying input and output 
/// files, aligners, classifiers, and various other parameters.
#[derive(Args, Debug)]
#[command(group(
    ArgGroup::new("aligner_classifier")
        .required(true)
        .args(&["aligner", "classifier"]), 
))]
#[command(group(
    ArgGroup::new("aligner_classifier_index")
        .required(true)
        .args(&["aligner_index", "classifier_index"]), 
))]
pub struct CleanArgs {
    /// Input read files (optional .gz)
    ///
    /// One or two input read files, can be in gzipped format. This parameter is required and multiple file
    /// can be specified (1 for long reads or 2 for paired-end short reads) either consecutively or using 
    /// multiple input arguments, for example: `-i R1.fq.gz -i R2.fq.gz` or `-i R1.fq.gz R2.fq.gz`
    #[arg(short, long, num_args(0..))]
    input: Vec<PathBuf>,
    /// Output read files (optional .gz)
    ///
    /// One or two output read files, can be in gzipped format. This parameter is required and multiple 
    /// files can be specified either consecutively or using multiple output arguments
    /// for example: `-o R1.fq.gz -o R2.fq.gz` or `-o R1.fq.gz R2.fq.gz`
    #[arg(short, long, num_args(0..))]
    output: Vec<PathBuf>,
    /// Read extraction instead of depletion
    ///
    /// Enable this option to extract reads matching the specified criteria instead 
    /// of depleting them.
    #[arg(short, long)]
    reverse: bool,
    /// Aligner to use
    ///
    /// Aligner to be used for the cleaning process. Options include 
    /// bowtie2, minimap2, and strobealign. If compiled with the `mm2`
    /// feature, the integrated minimap2-rs aligner becomes available.
    #[arg(long, short)]
    aligner: Option<Aligner>,
    /// Aligner index file
    ///
    /// Path to the aligner index file. This file is required for the 
    /// selected aligner. Can be a FASTA (optional .gz) for strobealign
    /// or minimap2 aligners.
    #[arg(long, short='A')]
    aligner_index: Option<PathBuf>,
    /// Classifier to use
    ///
    /// Classifier to be used for the cleaning process. Options include 
    /// Kraken2 and Metabuli.
    #[arg(long, short)]
    classifier: Option<Classifier>,
    /// Classifier index directory
    ///
    /// Path to the classifier index directory (Kraken2, Metabuli)
    #[arg(long, short='C')]
    classifier_index: Option<PathBuf>,
    /// Taxa and all sub-taxa to deplete using classifiers
    ///
    /// List of taxa names or taxids. All reads associated with these 
    /// taxa and their sub-taxa will be depleted.
    #[arg(long, short='T', num_args(0..))]
    taxa: Vec<String>,
    /// Taxa to deplete directly using classifiers
    ///
    /// List of taxa names or taxids to be directly depleted without 
    /// considering sub-taxa.
    #[arg(long, short='D', num_args(0..))]
    taxa_direct: Vec<String>,
    /// Number of threads to use for aligner and classifier
    ///
    /// Number of threads to be used by the aligner and classifier. 
    /// The default value is 4.
    #[arg(short, long, default_value = "4")]
    threads: usize,
    /// Summary output file (.json)
    ///
    /// Path to a JSON file for storing summary information about the 
    /// cleaning process.
    #[arg(short, long)]
    json: Option<PathBuf>,
    /// Optional working directory
    ///
    /// Working directory for temporary files. If not provided, the system 
    /// temporary directory will be used.
    #[arg(short, long)]
    workdir: Option<PathBuf>,
    /// Read identifier file (.tsv)
    ///
    /// Path to a TSV file containing read identifiers. This file can 
    /// be used to identify reads that were depleted or extracted.
    #[arg(long, short='R')]
    read_ids: Option<PathBuf>,
}
impl CleanArgs {
    /// Validates the provided arguments and builds a `Scrubby` instance.
    ///
    /// This method checks the provided arguments for consistency and constructs 
    /// a `Scrubby` instance based on the validated arguments.
    ///
    /// # Returns
    ///
    /// * `Result<Scrubby, ScrubbyError>` - Ok with the constructed Scrubby instance, otherwise an error.
    /// 
    /// # Example
    ///
    /// ```
    /// let clean_args = CleanArgs::parse();
    /// let scrubby = clean_args.validate_and_build().unwrap();
    /// ```
    pub fn validate_and_build(&self) -> Result<Scrubby, ScrubbyError> {
        let scrubby = ScrubbyBuilder::new(
            self.input.clone(), 
            self.output.clone()
        )
            .json(self.json.clone())
            .workdir(self.workdir.clone())
            .read_ids(self.read_ids.clone())
            .reverse(self.reverse)
            .threads(self.threads)
            .aligner(self.aligner.clone()) 
            .classifier(self.classifier.clone())
            .aligner_index(self.aligner_index.clone()) 
            .classifier_index(self.classifier_index.clone()) 
            .taxa(self.taxa.clone())
            .taxa_direct(self.taxa_direct.clone())
            .build()?;

        Ok(scrubby)
    }
}


#[derive(Args, Debug)]
pub struct ClassifierArgs {
    /// Input read files (can be compressed with .gz)
    ///
    /// One or two input read files. These files can be in gzipped format.
    /// This parameter is required and multiple files can be specified (1 for long
    /// reads or 2 for paired-end short reads) either consecutively or using multiple
    /// input arguments, for example: `-i R1.fq.gz -i R2.fq.gz` or `-i R1.fq.gz R2.fq.gz`
    #[arg(short, long, num_args(0..))]
    input: Vec<PathBuf>,
    /// Output read files (can be compressed with .gz)
    ///
    /// One or two output read files. These files will store the processed 
    /// data and can be in gzipped format. This parameter is required and multiple 
    /// files can be specified either consecutively or using multiple output arguments
    /// for example: `-o R1.fq.gz -o R2.fq.gz` or `-o R1.fq.gz R2.fq.gz`
    #[arg(short, long, num_args(0..))]
    output: Vec<PathBuf>,
    /// Read extraction instead of depletion
    ///
    /// Enable this option to extract reads matching the specified criteria instead 
    /// of depleting them.
    #[arg(short, long)]
    reverse: bool,
    /// Kraken-style report output from classifier
    ///
    /// Specify the path to the Kraken-style report file generated by the classifier.
    /// This file will contain the taxonomy results and summary metrics.
    #[arg(short='k', long)]
    report: PathBuf,
    /// Kraken-style read classification output
    ///
    /// Provide the path to the Kraken-style read classification file. This file 
    /// will contain the detailed read classifications from the classifier.
    #[arg(short='j', long)]
    reads: PathBuf,
    /// Classifier output style
    ///
    /// Classifier to be used for the cleaning process. Classifier that can  
    /// be used for this purpose are Kraken2 and Metabuli.
    #[arg(short, long)]
    classifier: Classifier,
    /// Taxa and all sub-taxa to deplete using classifiers
    ///
    /// List of taxa names or taxids. All reads associated with these 
    /// taxa and their sub-taxa will be depleted.
    #[arg(long, short='T', num_args(0..))]
    taxa: Vec<String>,
    /// Taxa to deplete directly using classifiers
    ///
    /// List of taxa names or taxids to be directly depleted without 
    /// considering sub-taxa.
    #[arg(long, short='D', num_args(0..))]
    taxa_direct: Vec<String>,
    /// Summary output file (.json)
    ///
    /// Path to a JSON file for storing summary information about the 
    /// cleaning process.
    #[arg(short, long)]
    json: Option<PathBuf>,
    /// Optional working directory
    ///
    /// Working directory for temporary files. If not provided, the system 
    /// temporary directory will be used.
    #[arg(short, long)]
    workdir: Option<PathBuf>,
    /// Read identifier file (.tsv)
    ///
    /// Path to a TSV file containing read identifiers. This file will 
    /// be used to identify specific reads for depletion or extraction.
    #[arg(short='R', long)]
    read_ids: Option<PathBuf>,
}
impl ClassifierArgs {
    /// Validates the provided arguments and builds a `Scrubby` instance.
    ///
    /// This method checks the provided arguments for consistency and constructs 
    /// a `Scrubby` instance based on the validated arguments.
    ///
    /// # Returns
    ///
    /// * `Result<Scrubby, ScrubbyError>` - Ok with the constructed Scrubby instance, otherwise an error.
    /// 
    /// # Example
    ///
    /// ```
    /// let clean_args = ClassifierArgs::parse();
    /// let scrubby = clean_args.validate_and_build().unwrap();
    /// ```
    pub fn validate_and_build(&self) -> Result<Scrubby, ScrubbyError> {
        let scrubby = ScrubbyBuilder::new(
            self.input.clone(), 
            self.output.clone()
        )
            .json(self.json.clone())
            .workdir(self.workdir.clone())
            .read_ids(self.reads.clone())
            .reverse(self.reverse)
            .classifier(self.classifier.clone())
            .classifier_reads(self.reads.clone())
            .classifier_report(self.report.clone())
            .taxa(self.taxa.clone())
            .taxa_direct(self.taxa_direct.clone())
            .build_classifier()?;

        Ok(scrubby)
    }
}


#[derive(Args, Debug)]
pub struct AlignmentArgs {
    /// Input read files (can be compressed with .gz)
    ///
    /// One or two input read files. These files can be in gzipped format.
    /// This parameter is required and multiple files can be specified (1 for long
    /// reads or 2 for paired-end short reads) either consecutively or using multiple
    /// input arguments, for example: `-i R1.fq.gz -i R2.fq.gz` or `-i R1.fq.gz R2.fq.gz`
    #[arg(short, long, num_args(0..))]
    input: Vec<PathBuf>,
    /// Output read files (can be compressed with .gz)
    ///
    /// One or two output read files. These files will store the processed 
    /// data and can be in gzipped format. This parameter is required and multiple 
    /// files can be specified either consecutively or using multiple output arguments
    /// for example: `-o R1.fq.gz -o R2.fq.gz` or `-o R1.fq.gz R2.fq.gz`
    #[arg(short, long, num_args(0..))]
    output: Vec<PathBuf>,
    /// Read extraction instead of depletion
    ///
    /// Enable this option to extract reads matching the specified criteria instead 
    /// of depleting them.
    #[arg(short, long)]
    reverse: bool,
    /// Alignment file in SAM/BAM/PAF/TXT format
    ///
    /// Specify the path to an alignment in SAM/BAM/PAF format (.sam, .bam, .cram, .paf),  
    /// or a read identifier file for any reads to deplete directly (.txt). PAF format
    /// alignments and rad identifier files can be compressed (.gz, .xz, .bz).
    #[arg(short, long)]
    alignment: PathBuf,
    /// Explicit alignment format
    /// 
    /// Otherwise format is determined from alignment extension, supported
    /// extensions are: 'paf', 'sam', 'bam', 'cram' or 'txt', with or without
    /// compression extensions 'paf.{gz, xz, bz}' or 'txt.{gz, xz, bz}'
    #[arg(short, long)]
    format: Option<AlignmentFormat>,
    /// Minimum query alignment length filter.
    #[arg(short='l', long, default_value = "0")]
    min_len: u64,
    /// Minimum query alignment coverage filter.
    #[arg(short='c', long, default_value = "0")]
    min_cov: f64,
    /// Minimum mapping quality filter.
    #[arg(short='q', long, default_value = "0")]
    min_mapq: u8,
    /// Summary output file (.json)
    ///
    /// Path to a JSON file for storing summary information about the 
    /// cleaning process.
    #[arg(short, long)]
    json: Option<PathBuf>,
    /// Optional working directory
    ///
    /// Working directory for temporary files. If not provided, the system 
    /// temporary directory will be used.
    #[arg(short, long)]
    workdir: Option<PathBuf>,
    /// Read identifier file (.tsv)
    ///
    /// Path to a TSV file containing read identifiers. This file will 
    /// be used to identify specific reads for depletion or extraction.
    #[arg(short='R', long)]
    read_ids: Option<PathBuf>,
}
impl AlignmentArgs {
    /// Validates the provided arguments and builds a `Scrubby` instance.
    ///
    /// This method checks the provided arguments for consistency and constructs 
    /// a `Scrubby` instance based on the validated arguments.
    ///
    /// # Returns
    ///
    /// * `Result<Scrubby, ScrubbyError>` - Ok with the constructed Scrubby instance, otherwise an error.
    /// 
    /// # Example
    ///
    /// ```
    /// let clean_args = ClassifierArgs::parse();
    /// let scrubby = clean_args.validate_and_build().unwrap();
    /// ```
    pub fn validate_and_build(&self) -> Result<Scrubby, ScrubbyError> {
        let scrubby = ScrubbyBuilder::new(
            self.input.clone(), 
            self.output.clone()
        )
            .json(self.json.clone())
            .workdir(self.workdir.clone())
            .read_ids(self.read_ids.clone())
            .reverse(self.reverse)
            .alignment(self.alignment.clone())
            .min_query_length(self.min_len)
            .min_query_coverage(self.min_cov)
            .min_mapq(self.min_mapq)
            .build_alignment()?;

        Ok(scrubby)
    }
}


#[derive(Args, Debug)]
pub struct DownloadArgs {
    /// Index name to download 
    /// 
    /// Default is for 'Bowtie2' unless `--aligner` or
    /// `--classfier` are set explicitly.
    #[arg(short, long, required=true, num_args(0..))]
    pub index: Vec<ScrubbyIndex>,
    /// Output directory for index download
    /// 
    /// Output directory will be created if it does not exist.
    #[arg(short, long)]
    pub outdir: PathBuf,
    /// Download index for one or more aligners 
    #[arg(short, long, num_args(0..))]
    pub aligner: Option<Vec<Aligner>>,
    /// Download index for one or more classifiers 
    #[arg(short, long, num_args(0..))]
    pub classfier: Option<Vec<Classifier>>,
    /// List available index names and exit
    #[arg(short, long)]
    pub list: bool,
    /// Download timeout in minutes - increase for large files and slow connections
    #[arg(short, long, default_value="60")]
    pub timeout: u64,
}
impl DownloadArgs {
    /// Validates the provided arguments and builds a `ScrubbyDownloader` instance.
    ///
    /// This method checks the provided arguments for consistency and constructs 
    /// a `ScrubbyDownloader` instance based on the validated arguments.
    ///
    /// # Returns
    ///
    /// * `Result<ScrubbyDownloader, ScrubbyError>` - Ok with the constructed ScrubbyDownloader instance, otherwise an error.
    /// 
    /// # Example
    ///
    /// ```
    /// let clean_args = ClassifierArgs::parse();
    /// let scrubby_dl = clean_args.validate_and_build().unwrap();
    /// ```
    pub fn validate_and_build(&self) -> Result<ScrubbyDownloader, ScrubbyError> {
        
        let downloader = ScrubbyDownloaderBuilder::new(
            self.outdir.clone(), self.index.clone()
        )
        .classifier(self.classfier.clone())
        .aligner(self.aligner.clone())
        .build()?;

        Ok(downloader)
    }
}

/// Configures the styles for the command-line interface.
pub fn get_styles() -> clap::builder::Styles {
    clap::builder::Styles::styled()
        .usage(
            anstyle::Style::new()
                .bold()
                .underline()
                .fg_color(Some(anstyle::Color::Ansi(anstyle::AnsiColor::Yellow))),
        )
        .header(
            anstyle::Style::new()
                .bold()
                .underline()
                .fg_color(Some(anstyle::Color::Ansi(anstyle::AnsiColor::Yellow))),
        )
        .literal(
            anstyle::Style::new()
                .bold()
                .fg_color(Some(anstyle::Color::Ansi(anstyle::AnsiColor::Green))),
        )
}