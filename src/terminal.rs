use std::path::PathBuf;
use clap::{ArgGroup, Args, Parser, Subcommand};

use crate::scrubby::{Aligner, Classifier, Scrubby, ScrubbyBuilder};
use crate::error::ScrubbyError;

/// Scrubby: background depletion for clinical metagenomic diagnostics
///
/// Scrubby provides fast and informed choices for taxonomic background read
/// depletion (host, microbial, etc. - extraction with --reverse) from paired-end
/// short read (Illumina) or long read (ONT) metagenomic sequencing appplications.
#[derive(Debug, Parser)]
#[command(author, version, about)]
#[command(styles=get_styles())]
#[command(arg_required_else_help(true))]
#[clap(name = "scrubby", version)]
pub struct App {
    /// Output logs to file instead of terminal
    ///
    /// Specify a file path to store the output logs. If not provided,
    /// logs will be displayed on the terminal.
    ///
    /// # Example
    ///
    /// ```
    /// scrubby --log-file logs.txt clean --input input.fastq --output output.fastq --aligner bowtie2 --aligner-index bowtie2_index
    /// ```
    #[arg(short, long)]
    pub log_file: Option<PathBuf>,

    #[clap(subcommand)]
    pub command: Commands,
}

/// Enumeration of available commands for Scrubby.
#[derive(Debug, Subcommand)]
pub enum Commands {
    /// Clean metagenomic data based on specified parameters.
    Clean(CleanArgs),
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
    /// Input read files (can be compressed with .gz)
    ///
    /// One or two input read files. These files can be in gzipped format.
    /// This parameter is required and multiple files can be specified (1 for long
    /// reads or 2 for paired-end short reads).
    #[arg(short, long, num_args(0..))]
    input: Vec<PathBuf>,
    /// Output read files (can be compressed with .gz)
    ///
    /// One or two output read files. These files will store the processed 
    /// data and can be in gzipped format. This parameter is required and multiple 
    /// files can be specified.
    #[arg(short, long, num_args(0..))]
    output: Vec<PathBuf>,
    /// Aligner to use
    ///
    /// Aligner to be used for the cleaning process. Options include 
    /// Bowtie2, Minimap2, and Strobealign.
    #[arg(long, short)]
    aligner: Option<Aligner>,
    /// Aligner index file
    ///
    /// Path to the aligner index file. This file is required for the 
    /// selected aligner.
    #[arg(long)]
    aligner_index: Option<PathBuf>,
    /// Classifier to use
    ///
    /// Classifier to be used for the cleaning process. Options include 
    /// Kraken2 and Metabuli.
    #[arg(long, short)]
    classifier: Option<Classifier>,
    /// Classifier index file
    ///
    /// Path to the classifier index file. This file is required for the 
    /// selected classifier.
    #[arg(long)]
    classifier_index: Option<PathBuf>,
    /// Taxa and all sub-taxa to deplete using classifiers
    ///
    /// List of taxa names or taxids. All reads associated with these 
    /// taxa and their sub-taxa will be depleted.
    #[arg(long, num_args(0..))]
    taxa: Vec<String>,
    /// Taxa to deplete directly using classifiers
    ///
    /// List of taxa names or taxids to be directly depleted without 
    /// considering sub-taxa.
    #[arg(long, num_args(0..))]
    taxa_direct: Vec<String>,
    /// Read identifier file (.tsv)
    ///
    /// Path to a TSV file containing read identifiers. This file will 
    /// be used to identify specific reads for depletion or extraction.
    #[arg(short, long)]
    reads: Option<PathBuf>,
    /// Summary output file (.json)
    ///
    /// Path to a JSON file for storing summary information about the 
    /// cleaning process.
    #[arg(short, long)]
    json: Option<PathBuf>,
    /// Optional working directory
    ///
    /// Wrking directory for temporary files. If not provided, the system 
    /// temporary directory will be used.
    #[arg(short, long)]
    workdir: Option<PathBuf>,
    /// Read extraction instead of depletion
    ///
    /// Enable this option to extract reads matching the specified criteria instead 
    /// of depleting them.
    #[arg(short, long)]
    reverse: bool,
    /// Number of threads to use for aligner and classifier
    ///
    /// Number of threads to be used by the aligner and classifier. 
    /// The default value is 4.
    #[arg(short, long, default_value = "4")]
    threads: usize,
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
        let scrubby = ScrubbyBuilder::new(self.input.clone(), self.output.clone())
                .reads(self.reads.clone())
                .json(self.json.clone())
                .workdir(self.workdir.clone())
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

/// Configures the styles for the command-line interface.
pub fn get_styles() -> clap::builder::Styles {
    clap::builder::Styles::styled()
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
