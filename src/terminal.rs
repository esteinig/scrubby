use std::path::PathBuf;
use clap::{ArgGroup, Args, Parser, Subcommand};

use crate::scrubby::{Aligner, Classifier, Scrubby, ScrubbyBuilder};
use crate::error::ScrubbyError;

/// Scrubby: background depletion for clinical metagenomic diagnostics
#[derive(Debug, Parser)]
#[command(author, version, about)]
#[command(styles=get_styles())]
#[command(arg_required_else_help(true))]
#[clap(name = "scrubby", version)]

// Global options and commands
pub struct App {
    #[clap(subcommand)]
    pub command: Commands,
    /// Output logs to file instead of terminal
    #[arg(short, long)]
    pub log_file: Option<PathBuf>,
}


/* Commands */

#[derive(Debug, Subcommand)]
pub enum Commands {
    Clean(CleanArgs),
}


/// Command-line arguments for the cleaning operation
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
    /// Input read files {gz}
    #[arg(short, long, num_args(0..))]
    input: Vec<PathBuf>,
    /// Output read files {gz}
    #[arg(short, long, num_args(0..))]
    output: Vec<PathBuf>,
    /// Aligner to use
    #[arg(long, short)]
    aligner: Option<Aligner>,
    /// Aligner index file
    #[arg(long)]
    aligner_index: Option<PathBuf>,
    /// Classifier to use
    #[arg(long, short)]
    classifier: Option<Classifier>,
    /// Classifier index file
    #[arg(long)]
    classifier_index: Option<PathBuf>,
    /// Taxa and all sub-taxa to deplete using classifiers
    #[arg(long, num_args(0..))]
    taxa: Vec<String>,
    /// Taxa to deplete directly using classifiers
    #[arg(long, num_args(0..))]
    taxa_direct: Vec<String>,
    /// Read identifier file (.tsv)
    #[arg(short, long)]
    reads: Option<PathBuf>,
    /// Summary output file (.json)
    #[arg(short, long)]
    json: Option<PathBuf>,
    /// Optional working directory
    #[arg(short, long)]
    workdir: Option<PathBuf>,
    /// Read extraction instead of depletion
    #[arg(short, long)]
    reverse: bool,
    /// Keep intermediate files
    #[arg(short, long)]
    keep: bool,
    /// Allow unpaired read depletion
    #[arg(short, long)]
    unpaired: bool,
    /// Number of threads to use
    #[arg(short, long, default_value = "4")]
    threads: usize,
}

impl CleanArgs {
    pub fn validate_and_build(&self) -> Result<Scrubby, ScrubbyError> {

        let scrubby = ScrubbyBuilder::new(self.input.clone(), self.output.clone())
                .reads(self.reads.clone())
                .json(self.json.clone())
                .workdir(self.workdir.clone())
                .reverse(self.reverse)
                .keep(self.keep)
                .unpaired(self.unpaired)
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
