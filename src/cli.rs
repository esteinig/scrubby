use std::ffi::{OsString, OsStr};
use std::path::PathBuf;
use structopt::StructOpt;
use structopt::clap::AppSettings;
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

/// MGP-DEPLETE command-line interface
#[derive(Debug, StructOpt)]
pub struct Cli {
    #[structopt(subcommand)]
    pub commands: Commands,
}

#[derive(Debug, StructOpt)]
pub enum Commands {
    #[structopt(global_settings = &[AppSettings::ColoredHelp, AppSettings::ArgRequiredElseHelp])]
    /// Clean seqeunce data by removing background taxa (k-mer) or host reads (alignment)
    Scrub {
        /// Input filepath(s) (fa, fq, gz, bz).
        ///
        /// For paired Illumina you may either pass this flag twice `-i r1.fq -i r2.fq` or give two
        /// files consecutively `-i r1.fq r2.fq`. Read identifiers for paired-end Illumina reads
        /// are assumed to be the same in forward and reverse read files (modern format) without trailing
        /// read orientations `/1` or `/2`.
        #[structopt(
            short,
            long,
            parse(try_from_os_str = check_file_exists),
            multiple = true,
            required = true
        )]
        input: Vec<PathBuf>,
        /// Output filepath(s) with contaminated reads removed.
        ///
        /// For paired Illumina you may either pass this flag twice `-o r1.fq -o r2.fq` or give two
        /// files consecutively `-o r1.fq r2.fq`. NOTE: The order of the pairs is assumed to be the
        /// same as that given for --input.
        #[structopt(
            short, 
            long, 
            parse(from_os_str), 
            multiple = true, 
            required = true
        )]
        output: Vec<PathBuf>,
        /// Extract reads instead of removing them (--output)
        ///
        /// This flagreverses the depletion and makes the command an extraction process 
        /// of reads that would otherwise be removed during depletion.
        #[structopt(short, long)]
        extract: bool,
        /// Kraken2 database directory path(s).
        ///
        /// Specify the path to the Kraken2 database directory. This only needs to be specified if you would like to
        /// run the Kraken2 analysis; otherwise `--kraken-report` and `--kraken-read` can be used. Note that multiple
        /// databases can be specified with `--kraken-db` which will be run in the order in which they are provided.
        /// You may either pass this flag twice `-k db1/ -k db2.` or give two files consecutively `-k db1/ db2/`.
        /// If multiple databases are provided, their names must be unique.
        /// 
        #[structopt(short = "k", long, parse(try_from_os_str = check_file_exists), multiple = true, required = true)]
        kraken_db: Vec<PathBuf>,
        /// Threads to use for Kraken2
        ///
        /// Specify the number of threads with which to run Kraken2.
        #[structopt(short = "j", long, default_value = "4")]
        kraken_threads: u32,
        /// Taxa to deplete from reads classified with Kraken2.
        ///
        /// You may specify multiple taxon names or taxonomic identifiers by passing this flag
        /// multiple times `-t Archaea -t 9606` or give taxa consecutively `-t Archaea 9606`.
        /// Kraken reports are parsed and every taxonomic level below the provided taxa will
        /// be provided for depletion of reads classified at these levels. For example, when
        /// providing `Archaea` (Domain) all taxonomic levels below the Domain level are removed
        /// until the next level of the same rank or higher is encountere in the report. This means
        /// that generally, higher levels than Domain should be specified with `--kraken-taxa-direct`.
        /// For example, if the level `cellular organisms` (R2) should be specified with `--kraken-taxa`
        /// it may also remove `Viruses` (D) if it is located below the `cellular organisms` level in
        /// the report. 
        #[structopt(short = "t", long, multiple = true, required = false)]
        kraken_taxa: Vec<String>,
        /// Taxa to deplete directly from reads classified with Kraken2.
        ///
        /// Additional taxon names or taxonomic identifiers at levels above the domain level provided with
        /// `--kraken-taxa` can be specified with this argument. These are directly to the list of taxons
        /// and not parsed from the report. For example, to retain Viruses one can specify the domains
        /// `-t Archaea -t Bacteria -t Eukaryota` with `--kraken-taxa` and add 
        /// `-d 'other sequences' -d 'cellular organsisms' -d root` with `--kraken-taxa-direct`.
        #[structopt(short = "d", long, multiple = true, required = false)]
        kraken_taxa_direct: Vec<String>,
        /// Working directory containing intermediary files.
        /// 
        /// Path to a working directory which contains the alignment and intermediary output files
        /// from the programs called during scrubbing.
        #[structopt(short, long, parse(from_os_str))]
        workdir: Option<PathBuf>,
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
        /// Compression level to use if compressing output.
        #[structopt(
            short = "l",
            long,
            parse(try_from_str = parse_level),
            default_value="6",
            value_name = "1-9"
        )]
        compression_level: niffler::Level,
    }
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
            Commands::Scrub { input, output, .. } => {
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

