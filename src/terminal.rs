use std::path::PathBuf;
use clap::{crate_version, Args, Parser, Subcommand};

use crate::prelude::*;
use crate::error::ScrubbyError;
use crate::utils::{ReadDifference, ReadDifferenceBuilder};

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
    /// Deplete or extract reads using aligners or classifiers.
    Reads(ReadsArgs),
    /// Deplete or extract reads from classifier outputs (Kraken2, Metabuli).
    Classifier(ClassifierArgs),
    /// Deplete or extract reads from aligner output with additional filters (SAM/BAM/PAF/GAF).
    Alignment(AlignmentArgs),
    /// List available indices and download files for aligners and classfiers.
    Download(DownloadArgs),
    /// Get read counts and identifiers of the difference between input and output read files.
    Diff(DiffArgs),
    /// Train and test the neural network for identity prediction.
    Nn(NeuralNetArgs)
}

/// Command-line arguments for the cleaning operation
///
/// This struct defines the arguments required for performing the cleaning 
/// operation using Scrubby. It includes options for specifying input and output 
/// files, aligners, classifiers, and various other parameters.
#[derive(Args, Debug)]
pub struct ReadsArgs {
    /// Input read files (optional .gz)
    ///
    /// One or two input read files, can be in gzipped format. This parameter is required and multiple file
    /// can be specified (1 for long reads or 2 for paired-end short reads) either consecutively or using 
    /// multiple input arguments, for example: '-i R1.fq.gz -i R2.fq.gz' or '-i R1.fq.gz R2.fq.gz'
    #[arg(short, long, num_args(0..))]
    input: Vec<PathBuf>,
    /// Output read files (optional .gz)
    ///
    /// One or two output read files, can be in gzipped format. This parameter is required and multiple 
    /// files can be specified either consecutively or using multiple output arguments for example:
    /// '-o R1.fq.gz -o R2.fq.gz' or '-o R1.fq.gz R2.fq.gz'. Output must be directed to files if 
    /// '--json' or '--read-ids' arguments are provided.
    #[arg(short, long, num_args(0..))]
    output: Vec<PathBuf>,
    /// Reference index for aligner or classifier
    ///
    /// Depending on whether --aligner or --classifier is chosen, the index is an 
    /// alignment index for 'bowtie2' (index), 'minimap2' and 'strobealign' 
    /// (index or FASTA) and 'minigraph' (graph index or FASTA) or a classifier
    /// index directory for Kraken2 (index) and Metabuli (index).
    #[arg(long, short='I')]
    index: PathBuf,
    /// Aligner to use, default is 'bowtie2' (paired) or 'minimap2' (single)
    ///
    /// Aligner to be used for the cleaning process. Default for paired-end short 
    /// reads is 'bowtie2'. Default for long reads is 'minimap2'. If compiled with 
    /// 'mm2' feature, the integrated 'minimap2-rs' aligner becomes available and
    /// is the default for both short and long reads.
    #[arg(long, short)]
    aligner: Option<Aligner>,
    /// Minimap2 or minigraph preset to use, defaults are: 'sr' | 'map-ont' | 'lr'
    ///
    /// Default is 'sr' for paired-end short reads, 'map-ont' for long reads
    /// and 'lr' for long reads with 'minigraph'.
    #[arg(long, short)]
    preset: Option<Preset>,
    /// Classifier to use
    ///
    /// Classifier to be used for the cleaning process. Options include 
    /// Kraken2 and Metabuli. Classifiers
    #[arg(long, short)]
    classifier: Option<Classifier>,
    /// Taxa and all sub-taxa to deplete using classifiers
    ///
    /// List of taxon names (case sensitive) or taxonomic identifiers. All reads associated with these 
    /// taxa and their descendent taxa will be depleted. Requires the taxon to be present in the index
    /// taxonomy. For example '--taxa Staphylococcus' would deplete all reads classified at genus level
    /// which includes all reads classfied at species level and below.
    #[arg(long, short='T', num_args(0..))]
    taxa: Vec<String>,
    /// Taxa to deplete directly using classifiers
    ///
    /// List of taxa names or taxonomic identifiers to be directly depleted without 
    /// considering descendent taxa. For example '--taxa-direct 9606' would deplete
    /// only reads directly classified as 'Homo sapiens' at species level.
    #[arg(long, short='D', num_args(0..))]
    taxa_direct: Vec<String>,
    /// Additional aligner arguments
    ///
    /// Aligner arguments must be a quoted string e.g. '-m 40'
    #[arg(long, short='A', allow_hyphen_values=true)]
    aligner_args: Option<String>,
    /// Additional classifier arguments
    ///
    /// Classifier arguments must be a quoted string e.g. '--min-score 0.008'
    #[arg(long, short='C', allow_hyphen_values=true)]
    classifier_args: Option<String>,
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
    #[arg(short, long)]
    read_ids: Option<PathBuf>,
    /// Read extraction instead of depletion
    ///
    /// Enable this option to extract reads matching the specified criteria instead 
    /// of depleting them.
    #[arg(short, long)]
    extract: bool,
}
impl ReadsArgs {
    /// Validates the provided arguments and builds a 'Scrubby' instance.
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
    /// use clap::Parser;
    /// 
    /// let reads_Args = ReadsArgs::parse();
    /// let scrubby = reads_Args.validate_and_build().unwrap();
    /// ```
    pub fn validate_and_build(self) -> Result<Scrubby, ScrubbyError> {

        let command = std::env::args().collect::<Vec<String>>().join(" ");
        
        let builder = ScrubbyBuilder::new(
            self.input, 
            self.output
        )
            .command(command)
            .json(self.json)
            .workdir(self.workdir)
            .read_ids(self.read_ids)
            .extract(self.extract)
            .threads(self.threads)
            .index(self.index)
            .aligner(self.aligner)
            .classifier(self.classifier)
            .taxa(self.taxa)
            .taxa_direct(self.taxa_direct)
            .classifier_args(self.classifier_args)
            .aligner_args(self.aligner_args)
            .preset(self.preset);

        let scrubby = builder.build()?;

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
    /// input arguments, for example: '-i R1.fq.gz -i R2.fq.gz' or '-i R1.fq.gz R2.fq.gz'
    #[arg(short, long, num_args(0..))]
    input: Vec<PathBuf>,
    /// Output read files (can be compressed with .gz)
    ///
    /// One or two output read files. These files will store the processed 
    /// data and can be in gzipped format. This parameter is required and multiple 
    /// files can be specified either consecutively or using multiple output arguments
    /// for example: '-o R1.fq.gz -o R2.fq.gz' or '-o R1.fq.gz R2.fq.gz'. Output must be 
    /// directed to files if '--json' or '--read-ids' arguments are provided.
    #[arg(short, long, num_args(0..))]
    output: Vec<PathBuf>,
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
    #[arg(short, long)]
    read_ids: Option<PathBuf>,
    /// Read extraction instead of depletion
    ///
    /// Enable this option to extract reads matching the specified criteria instead 
    /// of depleting them.
    #[arg(short, long)]
    extract: bool,
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
    /// use clap::Parser;
    /// 
    /// let class_args = ClassifierArgs::parse();
    /// let scrubby = class_args.validate_and_build().unwrap();
    /// ```
    pub fn validate_and_build(self) -> Result<Scrubby, ScrubbyError> {

        let command = std::env::args().collect::<Vec<String>>().join(" ");

        let scrubby = ScrubbyBuilder::new(
            self.input, 
            self.output
        )
            .command(command)
            .json(self.json)
            .workdir(self.workdir)
            .read_ids(self.read_ids)
            .extract(self.extract)
            .classifier(self.classifier)
            .reads(self.reads)
            .report(self.report)
            .taxa(self.taxa)
            .taxa_direct(self.taxa_direct)
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
    /// input arguments, for example: '-i R1.fq.gz -i R2.fq.gz' or '-i R1.fq.gz R2.fq.gz'
    #[arg(short, long, num_args(0..))]
    input: Vec<PathBuf>,
    /// Output read files (can be compressed with .gz)
    ///
    /// One or two output read files. These files will store the processed 
    /// data and can be in gzipped format. This parameter is required and multiple 
    /// files can be specified either consecutively or using multiple output arguments
    /// for example: '-o R1.fq.gz -o R2.fq.gz' or '-o R1.fq.gz R2.fq.gz'. Output must be 
    /// to file if '--json' or '--read-ids' arguments are provided.
    #[arg(short, long, num_args(0..))]
    output: Vec<PathBuf>,
    /// Alignment file in PAF/GAF/TXT or SAM/BAM/CRAM, if compiled with 'htslib' feature.
    ///
    /// Specify the path to an alignment in SAM/BAM/CRAM/PAF/GAF format (.sam, .bam, .cram, .paf, .gaf),  
    /// or a read identifier file for any reads to deplete directly (.txt). PAF/GAF/TXT format
    /// can be compressed (.gz, .xz, .bz). Allows '-' to read from stdin, but input stream cannot be 
    /// compressed and requires explicit setting of '--format'.
    #[arg(short, long)]
    alignment: PathBuf,
    /// Explicit alignment format
    /// 
    /// Otherwise format is determined from alignment extension, supported
    /// extensions with or without compression extension are: 
    /// '{paf, gaf, txt}.{gz, xz, bz, bz2}' and '{sam, bam, cram}' if compiled
    /// with 'htslib' feature.
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
    #[arg(short, long)]
    read_ids: Option<PathBuf>,
    /// Read extraction instead of depletion
    ///
    /// Enable this option to extract reads matching the specified criteria instead 
    /// of depleting them.
    #[arg(short, long)]
    extract: bool,
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
    /// use clap::Parser;
    /// 
    /// let aln_args = AlignmentArgs::parse();
    /// let scrubby = aln_args.validate_and_build().unwrap();
    /// ```
    pub fn validate_and_build(self) -> Result<Scrubby, ScrubbyError> {

        let command = std::env::args().collect::<Vec<String>>().join(" ");

        let scrubby = ScrubbyBuilder::new(
            self.input, 
            self.output,
        )   
            .command(command)
            .json(self.json)
            .workdir(self.workdir)
            .read_ids(self.read_ids)
            .extract(self.extract)
            .alignment(self.alignment)
            .alignment_format(self.format)
            .min_query_length(self.min_len)
            .min_query_coverage(self.min_cov)
            .min_mapq(self.min_mapq)
            .build_alignment()?;

        Ok(scrubby)
    }
}


#[derive(Args, Debug, Clone)]
pub struct DownloadArgs {
    /// Index name to download 
    /// 
    /// Default is 'bowtie2' aligner unless '--aligner' or
    /// '--classfier' arguments are set explicitly.
    #[arg(short, long, num_args(0..))]
    pub name: Vec<ScrubbyIndex>,
    /// Output directory for index download
    /// 
    /// Output directory will be created if it does not exist.
    #[arg(short, long, default_value=".")]
    pub outdir: PathBuf,
    /// Download index for one or more aligners 
    #[arg(short, long, num_args(0..))]
    pub aligner: Option<Vec<Aligner>>,
    /// Download index for one or more classifiers 
    #[arg(short, long, num_args(0..))]
    pub classfier: Option<Vec<Classifier>>,
    /// Download timeout in minutes - increase for large files and slow connections
    #[arg(short, long, default_value="360")]
    pub timeout: u64,
    /// List available index names and exit
    #[arg(short, long)]
    pub list: bool,
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
    /// use clap::Parser;
    /// 
    /// let dl_args = Download::parse();
    /// let dl = dl_args.validate_and_build().unwrap();
    /// ```
    pub fn validate_and_build(self) -> Result<ScrubbyDownloader, ScrubbyError> {
        
        let downloader = ScrubbyDownloaderBuilder::new(
            self.outdir, self.name
        )
        .classifier(self.classfier)
        .aligner(self.aligner)
        .timeout(self.timeout)
        .build()?;

        Ok(downloader)
    }
}


#[derive(Args, Debug)]
pub struct DiffArgs {
    /// Input read files (.gz | .xz | .bz)
    ///
    /// One or two input read files. These files can be in gzipped format.
    /// This parameter is required and multiple files can be specified (1 for long
    /// reads or 2 for paired-end short reads) either consecutively or using multiple
    /// input arguments, for example: '-i R1.fq.gz -i R2.fq.gz' or '-i R1.fq.gz R2.fq.gz'
    #[arg(short, long, num_args(0..))]
    input: Vec<PathBuf>,
    /// Output read files (.gz | .xz | .bz)
    ///
    /// One or two output read files. These files store the processed  data
    /// and can be in gzipped format. This parameter is required and multiple 
    /// files can be specified either consecutively or using multiple output 
    /// arguments for example: '-o R1.fq.gz -o R2.fq.gz' or '-o R1.fq.gz R2.fq.gz'
    #[arg(short, long, num_args(0..))]
    output: Vec<PathBuf>,
    /// Summary output file (.json)
    ///
    /// Path to a JSON file for storing summary information about the 
    /// difference between input and output files.
    #[arg(short, long)]
    json: Option<PathBuf>,
    /// Read identifier file (.tsv)
    ///
    /// Path to a TSV file containing read identifiers. This file will 
    /// be used to identify specific reads that constitute the difference
    /// between input and output files.
    #[arg(short, long)]
    read_ids: Option<PathBuf>,
}
impl DiffArgs {
    /// Validates the provided arguments and builds a `ReadDifference` instance.
    ///
    /// This method checks the provided arguments for consistency and constructs 
    /// a `ReadDifference` instance based on the validated arguments.
    ///
    /// # Returns
    ///
    /// * `Result<ReadDifference, ScrubbyError>` - Ok with the constructed ReadDifference instance, otherwise an error.
    /// 
    /// # Example
    ///
    /// ```
    /// use clap::Parser;
    /// 
    /// let diff_args = DifferenceArgs::parse();
    /// let diff = diff_args.validate_and_build().unwrap();
    /// ```
    pub fn validate_and_build(self) -> Result<ReadDifference, ScrubbyError> {
        
        Ok(ReadDifferenceBuilder::new(
            &self.input, 
            &self.output
        )
            .json(self.json)
            .read_ids(self.read_ids)
            .build()?)
    }
}



#[derive(Args, Debug)]
pub struct NeuralNetArgs {
    /// Input reads
    #[arg(short, long, num_args(0..))]
    pub fastq: Vec<PathBuf>,
    /// Model weights
    #[arg(short, long)]
    pub model_weights: PathBuf,
    /// Alignment
    #[arg(short, long)]
    pub alignment: Option<PathBuf>,
    /// Predict input reads with model
    #[arg(short, long)]
    pub predict: bool,
    /// Check GPU connect
    #[arg(short, long)]
    pub check: bool,
    /// Train model from input reads 
    #[arg(short, long)]
    pub train: bool,
    /// Epochs to train
    #[arg(short, long, default_value="10")]
    pub epochs: u64,
    /// Train with batch size
    #[arg(short, long, default_value="32")]
    pub batch_size: usize,
    /// CUDA device to use
    #[arg(short, long, default_value="0")]
    pub device: usize,
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
