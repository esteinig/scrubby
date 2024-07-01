use std::ffi::{OsStr, OsString};
use std::path::PathBuf;
use structopt::clap::AppSettings;
use structopt::StructOpt;
use thiserror::Error;

/// A collection of custom errors relating to the command line interface for this package.
#[derive(Error, Debug, PartialEq)]
pub enum CliError {
    /// Indicates that a string cannot be parsed into a [`CompressionFormat`](#compressionformat).
    #[error("{0} is not a valid output format")]
    CompressionFormat(String),
    /// Indicates that a string cannot be parsed into a [`CompressionLevel`](#compressionlevel).
    #[error("{0} is not a valid compression level (1-9)")]
    CompressionLevel(String),
    /// Indicates a bad combination of input and output files was passed.
    #[error("Bad combination of input and output files: {0}")]
    BadInputOutputCombination(String),
}

/// Scrubby command-line application
#[derive(Debug, StructOpt)]
pub struct Cli {
    #[structopt(short = "f", long)]
    pub force: bool,
    #[structopt(subcommand)]
    pub commands: Commands,
}

#[allow(clippy::enum_variant_names)]
#[derive(Debug, StructOpt)]
pub enum Commands {
    #[structopt(global_settings = &[AppSettings::ColoredHelp, AppSettings::ArgRequiredElseHelp])]
    /// Deplete or extract reads using k-mer classification (Kraken2) and/or alignments (Minimap2, Bowtie2, Strobealign)
    ScrubReads {
        /// Input filepath(s) (fa, fq, gz, bz).
        ///
        /// For paired Illumina you may either pass this flag twice `-i r1.fq -i r2.fq` or give two
        /// files consecutively `-i r1.fq r2.fq`. Read identifiers for paired-end Illumina reads
        /// are assumed to be the same in forward and reverse read files (modern format) without trailing
        /// read orientations `/1` or `/2`.
        #[structopt(
            short = "i",
            long,
            parse(try_from_os_str = check_file_exists),
            multiple = true,
            required = true
        )]
        input: Vec<PathBuf>,
        /// Output filepath(s) with reads removed or extracted.
        ///
        /// For paired Illumina you may either pass this flag twice `-o r1.fq -o r2.fq` or give two
        /// files consecutively `-o r1.fq r2.fq`. NOTE: The order of the pairs is assumed to be the
        /// same as that given for --input.
        #[structopt(
            short = "o",
            long,
            parse(from_os_str),
            multiple = true,
            required = true
        )]
        output: Vec<PathBuf>,
        /// Extract reads instead of removing them.
        ///
        /// This flag reverses the depletion and makes the command an extraction process
        /// of reads that would otherwise be removed during depletion.
        #[structopt(short = "e", long)]
        extract: bool,
        /// Kraken2 database directory path(s).
        ///
        /// Specify the path to the database directory to be used for classification with `Kraken2`. This only needs to be specified
        /// if you would like to run the `Kraken2` analysis; otherwise `--kraken-report` and `--kraken-read` can be used.
        /// Note that multiple databases can be specified with `--kraken-db` which will be run and reads depleted/extracted
        /// in the order with which the database files were provided. You may either pass this flag twice `-k db1/ -k db2/`
        /// or give two files consecutively `-k db1/ db2/`.
        #[structopt(short = "k", long, parse(try_from_os_str = check_file_exists), multiple = true, required = false)]
        kraken_db: Vec<PathBuf>,
        /// Threads to use for Kraken2.
        ///
        /// Specify the number of threads with which to run `Kraken2`.
        #[structopt(short = "j", long, default_value = "4")]
        kraken_threads: u32,
        /// Taxa and sub-taxa (Domain and below) to include.
        ///
        /// You may specify multiple taxon names or taxonomic identifiers by passing this flag
        /// multiple times `-t Archaea -t 9606` or give taxa consecutively `-t Archaea 9606`.
        /// `Kraken2` reports are parsed and every taxonomic level below the provided taxon level will
        /// be included. Only taxa or sub-taxa that have reads directly assigned to them will be parsed.
        /// For example, when providing `Archaea` (Domain) all taxonomic levels below the `Domain` level are
        /// included until the next level of the same rank or higher is encountered in the report. This means
        /// that higher levels than `Domain` should be specified with `--kraken-taxa-direct`.
        #[structopt(short = "t", long, multiple = true, required = false)]
        kraken_taxa: Vec<String>,
        /// Taxa to include directly from reads classified.
        ///
        /// Additional taxon names or taxonomic identifiers can be specified with this argument,
        /// such as those above the `Domain` level. These are directly added to the list of taxa to include
        /// while parsing the report without considering sub-taxa. For example, to retain `Viruses` one can
        /// specify the domains `-t Archaea -t Bacteria -t Eukaryota` with `--kraken-taxa` and add
        /// `-d 'other sequences' -d 'cellular organsisms' -d root` with `--kraken-taxa-direct`.
        #[structopt(short = "d", long, multiple = true, required = false)]
        kraken_taxa_direct: Vec<String>,
        /// Metabuli database directory path(s).
        ///
        /// Specify the path to the database directory to be used for classification with `Metabuli`. This only needs to be specified
        /// if you would like to run the `Metabuli` analysis; otherwise `--metabuli-report` and `--metabuli-read` can be used.
        /// Note that multiple databases can be specified with `--metabuli-db` which will be run and reads depleted/extracted
        /// in the order with which the database files were provided. You may either pass this flag twice `-m db1/ -m db2/`
        /// or give two files consecutively `-k db1/ db2/`.
        #[structopt(long, parse(try_from_os_str = check_file_exists), multiple = true, required = false)]
        metabuli_db: Vec<PathBuf>,
        /// Threads to use for Kraken2.
        ///
        /// Specify the number of threads with which to run `Kraken2`.
        #[structopt(long, default_value = "4")]
        metabuli_threads: u32,
        /// Reference sequence or index file(s) for `minimap2`.
        ///
        /// Specify the index file (.mmi) or the reference sequence(s) (.fasta) for alignment with `minimap2`. Note that creating
        /// the index file may take some time with larger genomes. Multiple references can be specified with `--minimap2-index` which
        /// will be run and reads depleted/extracted in the order with which the database files were provided. You may either pass this
        /// flag twice `-m idx1.mmi -m idx2.mmi` or give two files consecutively `-m idx1.mmi idx2.mmi`.
        #[structopt(short = "m", long, parse(try_from_os_str = check_file_exists), multiple = true, required = false)]
        minimap2_index: Vec<PathBuf>,
        /// Minimap2 preset configuration (sr|map-ont|map-hifi|map-pb).
        ///
        /// Specify the preset configuration for `minimap2`.
        #[structopt(
            short = "x",
            long,
            default_value = "sr",
            multiple = false,
            required = false,
            value_name = "sr|map-ont|map-hifi|map-pb",
            case_insensitive = true,
            hide_possible_values = true,
            possible_values = &["sr", "map-ont", "map-hifi", "map-pb"],
        )]
        minimap2_preset: String,
        /// Threads to use for `minimap2`.
        ///
        /// Specify the number of threads with which to run `minimap2`.
        #[structopt(short = "n", long, default_value = "4")]
        minimap2_threads: u32,
        /// Reference sequence or index file(s) for `strobealign`.
        ///
        /// Specify the index file (.sti) or the reference sequence(s) (.fasta) for alignment with `strobealign`. When a reference
        /// sequence is specified, the format is checked for being FASTA and a new index is created optimised for the input read length
        /// for each run, this may take time with larger genomes. If a precomputed index is supplied, the index read length should match the
        /// input read length (see strobealign manual). Note that multiple references can be specified with `--strobealign-index` which will
        /// be run and reads depleted/extracted in the order with which the database files were provided. You may either pass this flag twice
        /// `-s idx1.sti -s idx2.sti` or give two files consecutively `-s idx1.sti idx2.sti`.
        #[structopt(long, parse(try_from_os_str = check_file_exists), multiple = true, required = false)]
        strobealign_index: Vec<PathBuf>,
        /// Strobalign alignment mode (map|align).
        ///
        /// Specify the alignment mode for `strobelalign`. Mapping (`map`) ignores quality scores and outputs PAF,
        /// alignment (`align`) uses quality scores and outputs SAM.
        #[structopt(
            long,
            default_value = "align",
            multiple = false,
            required = false,
            value_name = "map|align",
            case_insensitive = true,
            hide_possible_values = true,
            possible_values = &["map", "align"],
        )]
        strobealign_mode: String,
        /// Threads to use for `minimap2`.
        ///
        /// Specify the number of threads with which to run `minimap2`.
        #[structopt(long, default_value = "4")]
        strobealign_threads: u32,
        /// Minimum query alignment length filter.
        #[structopt(short = "l", long, default_value = "0")]
        min_len: u64,
        /// Minimum query alignment coverage filter.
        #[structopt(short = "c", long, default_value = "0")]
        min_cov: f64,
        /// Minimum mapping quality filter.
        #[structopt(short = "q", long, default_value = "0")]
        min_mapq: u8,
        /// Working directory containing intermediary files.
        ///
        /// Path to a working directory which contains the alignment and intermediary output files
        /// from the programs called during scrubbing. By default is the working output directory
        /// is named with a timestamp in the format: `Scrubby_{YYYYMMDDTHHMMSS}`.
        #[structopt(short = "W", long, parse(from_os_str))]
        workdir: Option<PathBuf>,
        /// Output filepath for summary of depletion/extraction.
        ///
        /// This specified a JSON formatted output file that contains a summary of the
        /// depletion/extraction steps (number of reads, total/depleted/extracted/retained)
        #[structopt(short = "J", long, parse(from_os_str))]
        json: Option<PathBuf>,
        /// Keep the working directory and intermediate files.
        ///
        /// This flag specifies that we want to keep the working directory and all intermediate files;
        /// otherwise the working directory is deleted.
        #[structopt(short = "K", long)]
        keep: bool,
        /// Include read identifiers as a summary output file
        ///
        /// This flag specifies the output of a table (.tsv) with three columns and headers (id, ref, method) with 
        /// the read identifiers which were removed/extracted and the reference/method used to remove/extract them
        #[structopt(short = "R", long, parse(from_os_str))]
        reads: Option<PathBuf>,
        /// u: uncompressed; b: Bzip2; g: Gzip; l: Lzma
        ///
        /// Default is to attempt to infer the output compression format automatically from the filename
        /// extension (gz|bz|bz2|lzma). This option is used to override that.
        #[structopt(
            short = "O",
            long,
            value_name = "u|b|g|l",
            parse(try_from_str = parse_compression_format),
            possible_values = &["u", "b", "g", "l"],
            case_insensitive=true,
            hide_possible_values = true
        )]
        output_format: Option<niffler::compression::Format>,
        /// Compression level to use.
        #[structopt(
            short = "L",
            long,
            parse(try_from_str = parse_level),
            default_value="6",
            value_name = "1-9"
        )]
        compression_level: niffler::Level,
    },
    /// Deplete or extract reads using outputs from Kraken2
    ScrubKraken {
        /// Input filepath(s) (fa, fq, gz, bz).
        ///
        /// For paired Illumina you may either pass this flag twice `-i r1.fq -i r2.fq` or give two
        /// files consecutively `-i r1.fq r2.fq`. Read identifiers for paired-end Illumina reads
        /// are assumed to be the same in forward and reverse read files (modern format) without trailing
        /// read orientations `/1` or `/2`.
        #[structopt(
            short = "i",
            long,
            parse(try_from_os_str = check_file_exists),
            multiple = true,
            required = true
        )]
        input: Vec<PathBuf>,
        /// Output filepath(s) with reads removed or extracted.
        ///
        /// For paired Illumina you may either pass this flag twice `-o r1.fq -o r2.fq` or give two
        /// files consecutively `-o r1.fq r2.fq`. NOTE: The order of the pairs is assumed to be the
        /// same as that given for --input.
        #[structopt(
            short = "o",
            long,
            parse(from_os_str),
            multiple = true,
            required = true
        )]
        output: Vec<PathBuf>,
        /// Extract reads instead of removing them.
        ///
        /// This flag reverses the depletion and makes the command an extraction process
        /// of reads that would otherwise be removed during depletion.
        #[structopt(short = "e", long)]
        extract: bool,
        /// Kraken2 classified reads output.
        ///
        #[structopt(short = "k", long,  parse(try_from_os_str = check_file_exists), multiple = false, required = true)]
        kraken_reads: PathBuf,
        /// Kraken2 taxonomic report output.
        ///
        #[structopt(short = "r", long,  parse(try_from_os_str = check_file_exists), multiple = false, required = true)]
        kraken_report: PathBuf,
        /// Taxa and sub-taxa (Domain and below) to include.
        ///
        /// You may specify multiple taxon names or taxonomic identifiers by passing this flag
        /// multiple times `-t Archaea -t 9606` or give taxa consecutively `-t Archaea 9606`.
        /// `Kraken2` reports are parsed and every taxonomic level below the provided taxon level will
        /// be included. Only taxa or sub-taxa that have reads directly assigned to them will be parsed.
        /// For example, when providing `Archaea` (Domain) all taxonomic levels below the `Domain` level are
        /// included until the next level of the same rank or higher is encountered in the report. This means
        /// that higher levels than `Domain` should be specified with `--kraken-taxa-direct`.
        #[structopt(short = "t", long, multiple = true, required = false)]
        kraken_taxa: Vec<String>,
        /// Taxa to include directly from reads classified.
        ///
        /// Additional taxon names or taxonomic identifiers can be specified with this argument,
        /// such as those above the `Domain` level. These are directly added to the list of taxa to include
        /// while parsing the report without considering sub-taxa. For example, to retain `Viruses` one can
        /// specify the domains `-t Archaea -t Bacteria -t Eukaryota` with `--kraken-taxa` and add
        /// `-d 'other sequences' -d 'cellular organsisms' -d root` with `--kraken-taxa-direct`.
        #[structopt(short = "d", long, multiple = true, required = false)]
        kraken_taxa_direct: Vec<String>,
        /// Database name for JSON summary, by default uses --kraken-reads filestem
        ///
        /// This option provides an alternative name for the database in the JSON summary
        /// in cases where the input classification file is named e.g. {sample_id}.kraken
        /// which would not be particularly informativ in the summaries
        #[structopt(short = "n", long)]
        kraken_name: Option<String>,
        /// Working directory for intermediary files.
        ///
        /// Path to a working directory which contains the alignment and intermediary output files
        /// from the programs called during scrubbing. By default is the working output directory
        /// is named with a timestamp in the format: `Scrubby_{YYYYMMDDTHHMMSS}`.
        #[structopt(short = "W", long, parse(from_os_str))]
        workdir: Option<PathBuf>,
        /// Output filepath for summary of depletion/extraction.
        ///
        /// This specified a JSON formatted output file that contains a summary of the
        /// depletion/extraction steps (number of reads, total/depleted/extracted/retained)
        #[structopt(short = "J", long, parse(from_os_str))]
        json: Option<PathBuf>,
        /// u: uncompressed; b: Bzip2; g: Gzip; l: Lzma
        ///
        /// Default is to attempt to infer the output compression format automatically from the filename
        /// extension (gz|bz|bz2|lzma). This option is used to override that.
        #[structopt(
            short = "O",
            long,
            value_name = "u|b|g|l",
            parse(try_from_str = parse_compression_format),
            possible_values = &["u", "b", "g", "l"],
            case_insensitive=true,
            hide_possible_values = true
        )]
        output_format: Option<niffler::compression::Format>,
        /// Compression level to use.
        #[structopt(
            short = "L",
            long,
            parse(try_from_str = parse_level),
            default_value="6",
            value_name = "1-9"
        )]
        compression_level: niffler::Level,
    },
    /// Deplete or extract reads using alignments (PAF|SAM|BAM|CRAM)
    ScrubAlignment {
        /// Input filepath(s) (fa, fq, gz, bz).
        ///
        /// For paired Illumina you may either pass this flag twice `-i r1.fq -i r2.fq` or give two
        /// files consecutively `-i r1.fq r2.fq`. Read identifiers for paired-end Illumina reads
        /// are assumed to be the same in forward and reverse read files (modern format) without trailing
        /// read orientations `/1` or `/2`.
        #[structopt(
            short = "i",
            long,
            parse(try_from_os_str = check_file_exists),
            multiple = true,
            required = true
        )]
        input: Vec<PathBuf>,
        /// Output filepath(s) with reads removed or extracted.
        ///
        /// For paired Illumina you may either pass this flag twice `-o r1.fq -o r2.fq` or give two
        /// files consecutively `-o r1.fq r2.fq`. NOTE: The order of the pairs is assumed to be the
        /// same as that given for --input.
        #[structopt(
            short = "o",
            long,
            parse(from_os_str),
            multiple = true,
            required = true
        )]
        output: Vec<PathBuf>,
        /// Extract reads instead of removing them.
        ///
        /// This flag reverses the depletion and makes the command an extraction process
        /// of reads that would otherwise be removed during depletion.
        #[structopt(short = "e", long)]
        extract: bool,
        /// Alignment file (SAM/BAM/CRAM/PAF) or list of read identifiers (TXT)
        #[structopt(
            short = "a", long, parse(try_from_os_str = check_file_exists), required = true
        )]
        alignment: PathBuf,
        /// bam: SAM/BAM/CRAM alignment; paf: PAF alignment, txt: read identifiers
        ///
        /// Default is to attempt to infer the input alignment format automatically from the filename
        /// extension (.bam|.sam|.cram|.paf|.txt). This option is used to override that.
        #[structopt(
            short = "A",
            long,
            value_name = "bam|paf|txt|kraken",
            possible_values = &["bam", "paf", "txt"],
            case_insensitive=true,
            hide_possible_values=true
        )]
        alignment_format: Option<String>,
        /// Alignment name for JSON summary, by default uses --alignment filestem
        ///
        /// This option provides an alternative name for the alignment in the JSON summary
        /// in cases where the input alignment is named e.g. {sample_id}.paf which would
        /// not be particularly informativ in the summaries
        #[structopt(short = "n", long)]
        alignment_name: Option<String>,
        /// Minimum query alignment length filter.
        #[structopt(short = "l", long, default_value = "0")]
        min_len: u64,
        /// Minimum query alignment coverage filter.
        #[structopt(short = "c", long, default_value = "0")]
        min_cov: f64,
        /// Minimum mapping quality filter.
        #[structopt(short = "q", long, default_value = "0")]
        min_mapq: u8,
        #[structopt(short = "W", long, parse(from_os_str))]
        workdir: Option<PathBuf>,
        /// Output filepath for summary of depletion/extraction.
        ///
        /// This specified a JSON formatted output file that contains a summary of the
        /// depletion/extraction steps (number of reads, total/depleted/extracted/retained)
        #[structopt(short = "J", long, parse(from_os_str))]
        json: Option<PathBuf>,
        /// u: uncompressed; b: Bzip2; g: Gzip; l: Lzma
        ///
        /// Default is to attempt to infer the output compression format automatically from the filename
        /// extension (gz|bz|bz2|lzma). This option is used to override that.
        #[structopt(
            short = "O",
            long,
            value_name = "u|b|g|l",
            parse(try_from_str = parse_compression_format),
            possible_values = &["u", "b", "g", "l"],
            case_insensitive=true,
            hide_possible_values = true
        )]
        output_format: Option<niffler::compression::Format>,
        /// Compression level to use.
        #[structopt(
            short = "L",
            long,
            parse(try_from_str = parse_level),
            default_value="6",
            value_name = "1-9"
        )]
        compression_level: niffler::Level,
    },
}

impl Cli {
    /// Checks there is a valid and equal number of `--input` and `--output` arguments given.
    ///
    /// # Errors
    /// A [`CliError::BadInputOutputCombination`](#clierror) is returned for the following:
    /// - Either `--input` or `--output` are passed more than twice
    /// - An unequal number of `--input` and `--output` are passed
    pub fn validate_input_output_combination(&self) -> Result<(), CliError> {
        match &self.commands {
            Commands::ScrubReads { input, output, .. } => {
                let out_len = output.len();
                let in_len = input.len();
                if in_len > 2 {
                    let msg = String::from("Got more than 2 files for input.");
                    return Err(CliError::BadInputOutputCombination(msg));
                }
                if out_len > 2 {
                    let msg = String::from("Got more than 2 files for output.");
                    return Err(CliError::BadInputOutputCombination(msg));
                }
                if in_len != out_len {
                    let msg = format!("Got {} --input but {} --output", in_len, out_len);
                    return Err(CliError::BadInputOutputCombination(msg));
                }
            }
            Commands::ScrubKraken { .. } => {}
            Commands::ScrubAlignment { .. } => {}
        };
        Ok(())
    }
}

/// A utility function to validate whether an input files exist
fn check_file_exists(file: &OsStr) -> Result<PathBuf, OsString> {
    let path = PathBuf::from(file);
    let path_msg = format!("{:?} does not exist", path);
    if path.exists() {
        let abs_path = path.canonicalize().map_err(|_| OsString::from(path_msg))?;
        Ok(abs_path)
    } else {
        Err(OsString::from(path_msg))
    }
}

/// A utility function to validate compression format is in allowed values
fn parse_compression_format(s: &str) -> Result<niffler::compression::Format, CliError> {
    match s {
        "b" | "B" => Ok(niffler::Format::Bzip),
        "g" | "G" => Ok(niffler::Format::Gzip),
        "l" | "L" => Ok(niffler::Format::Lzma),
        "u" | "U" => Ok(niffler::Format::No),
        _ => Err(CliError::CompressionFormat(s.to_string())),
    }
}

/// A utility function to validate compression level is in allowed range
#[allow(clippy::redundant_clone)]
fn parse_level(s: &str) -> Result<niffler::Level, CliError> {
    let lvl = match s.parse::<u8>() {
        Ok(1) => niffler::Level::One,
        Ok(2) => niffler::Level::Two,
        Ok(3) => niffler::Level::Three,
        Ok(4) => niffler::Level::Four,
        Ok(5) => niffler::Level::Five,
        Ok(6) => niffler::Level::Six,
        Ok(7) => niffler::Level::Seven,
        Ok(8) => niffler::Level::Eight,
        Ok(9) => niffler::Level::Nine,
        _ => return Err(CliError::CompressionLevel(s.to_string())),
    };
    Ok(lvl)
}
