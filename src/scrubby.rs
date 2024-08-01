//! This module provides structures and implementations for Scrubby, 
//! which is used for read cleaning or extraction using k-mer classifiers
//! and aligners. It is primarily meant for benchmarking host depletion for
//! clinical metagenomic diagnostic applications in conjunction with syndromic
//! reference panel simulations and reference index construction including for
//! human and microbial pangenome databases with Cipher.
//! 
//! Scrubby provides optimised choices for clinical metagenomics short- and long
//! read data (Illumina, ONT) depending on target applications such as deep short read
//! sequencing of sterile sites (e.g. for central nervous system infections) where specificity 
//! is paramount, or on the other hand for targeted enrichment of data sensitive pathogens like
//! HIV sequencing from primary samples, where sensitivty and (near) complete removal of
//! host background for data storage or release becomes more important.
//! 
//! The primary structures and enumerations provided are `Aligner`, `Classifier`, `Preset`,
//! `Scrubby`, `ScrubbyConfig`, and `ScrubbyBuilder`. These components can be imported
//! through the prelude module: `scrubby::prelude::*`.

use serde::{Serialize, Deserialize};
use std::fs::create_dir_all;
use std::path::PathBuf;
use std::fmt;

use crate::cleaner::Cleaner;
use crate::error::ScrubbyError;
use crate::utils::IntoVecPathBuf;
use crate::report::ScrubbyReport;
use crate::alignment::AlignmentFormat;

/// Enum representing the available aligners.
#[derive(Serialize, Deserialize, Clone, Debug, PartialEq, clap::ValueEnum)]
pub enum Aligner {
    #[serde(rename="bowtie2")]
    Bowtie2,
    #[serde(rename="minimap2")]
    Minimap2,
    #[serde(rename="minigraph")]
    Minigraph,
    #[serde(rename="strobealign")]
    Strobealign,
    #[serde(rename="minimap2-rs")]
    #[cfg(feature = "mm2")]
    Minimap2Rs
}
impl Aligner {
    // Used for identification of pre-built-indices
    pub fn short_name(&self) -> &str {
        match self {
            Aligner::Bowtie2 => "bt",
            Aligner::Minimap2 => "mm",
            Aligner::Minigraph => "mm",
            Aligner::Strobealign => "st",
            #[cfg(feature = "mm2")]
            Aligner::Minimap2Rs => "mm"
        }
    }
}
impl fmt::Display for Aligner {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Aligner::Bowtie2 => write!(f, "bowtie2"),
            Aligner::Minimap2 => write!(f, "minimap2"),
            Aligner::Minigraph => write!(f, "minigraph"),
            Aligner::Strobealign => write!(f, "strobealign"),
            #[cfg(feature = "mm2")]
            Aligner::Minimap2Rs => write!(f, "minimap2-rs"),
        }
    }
}

/// Enum representing the available classifiers.
#[derive(Serialize, Deserialize, Clone, Debug, clap::ValueEnum)]
pub enum Classifier {
    #[serde(rename="kraken2")]
    Kraken2,
    #[serde(rename="metabuli")]
    Metabuli,
}

impl Classifier {
    pub fn short_name(&self) -> &str {
        match self {
            Classifier::Kraken2 => "k2",
            Classifier::Metabuli => "mb",
        }
    }
}
impl fmt::Display for Classifier {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Classifier::Kraken2 => write!(f, "kraken2"),
            Classifier::Metabuli => write!(f, "metabuli"),
        }
    }
}


/// TODO: Enum representing the available classifiers output styles
/// for direct classifier output cleaning 
#[derive(Serialize, Deserialize, Clone, Debug, clap::ValueEnum)]
pub enum ClassifierOutput {
    #[serde(rename="kraken2")]
    Kraken2,
    #[serde(rename="metabuli")]
    Metabuli,
    #[serde(rename="kraken2uniq")]
    Kraken2Uniq
}
impl fmt::Display for ClassifierOutput {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            ClassifierOutput::Kraken2 => write!(f, "kraken2"),
            ClassifierOutput::Metabuli => write!(f, "metabuli"),
            ClassifierOutput::Kraken2Uniq => write!(f, "kraken2uniq"),
        }
    }
}

#[derive(Serialize, Deserialize, Clone, Debug, PartialEq, clap::ValueEnum)]
pub enum Preset {
    LrHq,
    Splice,
    SpliceHq,
    Asm,        // minigraph + minimap
    Asm5,
    Asm10,
    Asm20,
    Sr,         // minigraph + minimap
    Lr,         // minigraph 
    MapPb,
    MapHifi,
    MapOnt,
    AvaPb,
    AvaOnt,
}
impl fmt::Display for Preset {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Preset::LrHq => write!(f, "lr:hq"),
            Preset::Splice => write!(f, "splice"),
            Preset::SpliceHq => write!(f, "splice:hq"),
            Preset::Asm => write!(f, "asm"),
            Preset::Asm5 => write!(f, "asm5"),
            Preset::Asm10 => write!(f, "asm10"),
            Preset::Asm20 => write!(f, "asm20"),
            Preset::Sr => write!(f, "sr"),
            Preset::Lr => write!(f, "lr"),
            Preset::MapPb => write!(f, "map-pb"),
            Preset::MapHifi => write!(f, "map-hifi"),
            Preset::MapOnt => write!(f, "map-ont"),
            Preset::AvaPb => write!(f, "ava-pb"),
            Preset::AvaOnt => write!(f, "ava-ont")
        }
    }
}

/// Main structure representing the Scrubby tool configuration.
#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct Scrubby {
    pub input: Vec<PathBuf>,
    pub output: Vec<PathBuf>,
    pub json: Option<PathBuf>,
    pub workdir: Option<PathBuf>,
    pub read_ids: Option<PathBuf>,
    pub extract: bool,
    pub keep: bool,
    pub threads: usize,
    pub config: ScrubbyConfig,
}

impl Scrubby {
    /// Creates a new `Scrubby` instance with default configuration. If both aligner and classifier are `None`,
    /// the default aligner (and preset) is selected depending on file input and output length: 
    /// 
    /// * 1 - long reads with `Aligner::Minimap2` and `Preset::MapOnt`
    /// * 2 - short reads with `Aligner::Bowtie2`
    /// 
    /// If compiled with the `mm2` feature, the following default aligners are selected:
    /// 
    /// * 1 - long reads with `Aligner::Minimap2` and `Preset::MapOnt`
    /// * 2 - short reads with `Aligner::Minimap2` and `Preset::Sr`
    ///
    /// # Arguments
    ///
    /// * `input` - A vector of input file paths.
    /// * `output` - A vector of output file paths.
    /// * `index` - An index file path.
    /// * `aligner` - Optional aligner configuration.
    /// * `classifier` - Optional classifier configuration.
    ///
    /// # Example
    ///
    /// ```
    /// use scrubby::{Scrubby, Aligner, Classifier};
    /// use std::path::PathBuf;
    ///
    /// let scrubby = Scrubby::new(
    ///     "input.fastq", 
    ///     "output.fastq", 
    ///     "kraken2_db/",
    ///     None, 
    ///     Classifier::Kraken2, 
    /// )).unwrap();
    /// ```
    pub fn new<P: IntoVecPathBuf, I: Into<PathBuf>, A: Into<Option<Aligner>>, C: Into<Option<Classifier>>>(
        input: P, 
        output: P, 
        index: I,
        aligner: A, 
        classifier: C, 
    ) -> Result<Self, ScrubbyError> {
        ScrubbyBuilder::new(input.into_vec_path_buf(), output.into_vec_path_buf())
            .index(index.into())
            .aligner(aligner)
            .classifier(classifier)
            .build()
    }

    /// Creates a new `ScrubbyBuilder` instance for constructing a `Scrubby` object. 
    ///
    /// # Arguments
    ///
    /// * `input` - A vector of input file paths.
    /// * `output` - A vector of output file paths.
    ///
    /// # Example
    ///
    /// ```
    /// use scrubby::Scrubby;
    /// use std::path::PathBuf;
    ///
    /// let builder = Scrubby::builder(
    ///     "input.fastq", 
    ///     "output.fastq"
    /// );
    /// ```
    pub fn builder<P: IntoVecPathBuf>(input: P, output: P) -> ScrubbyBuilder {
        ScrubbyBuilder::new(input, output)
    }

    /// Executes the cleaning process based on the Scrubby configuration.
    ///
    /// # Returns
    ///
    /// * `Result<(), ScrubbyError>` - Ok if the cleaning process completes successfully, otherwise an error.
    ///
    /// # Example
    ///
    /// ```
    /// use scrubby::Scrubby;
    /// 
    /// let scrubby = Scrubby::new(...).unwrap();
    /// scrubby.clean().unwrap();
    /// ```
    pub fn clean(&self) -> Result<(), ScrubbyError> {
        let cleaner = Cleaner::from_scrubby(self)?;

        if self.config.aligner.is_some() {
            cleaner.run_aligner()?;
        }
        if self.config.classifier.is_some() {
            cleaner.run_classifier()?;
        }
        if self.config.reads.is_some() && self.config.report.is_some() {
            cleaner.run_classifier_output()?;
        }
        if self.config.alignment.is_some() {
            cleaner.run_aligner_output()?;
        }
        if self.json.is_some() || self.read_ids.is_some() {
            ScrubbyReport::create(&self, true)?;
        }

        Ok(())
    }
}

/// Configuration structure for Scrubby
#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct ScrubbyConfig {
    pub aligner: Option<Aligner>,
    pub classifier: Option<Classifier>,
    pub index: Option<PathBuf>,
    pub aligner_index: Option<PathBuf>,
    pub alignment: Option<PathBuf>,
    pub classifier_index: Option<PathBuf>,
    pub reads: Option<PathBuf>,
    pub report: Option<PathBuf>,
    pub taxa: Vec<String>,
    pub taxa_direct: Vec<String>,
    pub classifier_args: Option<String>,
    pub aligner_args: Option<String>,
    pub unpaired: bool,
    pub paired_end: bool,
    pub samtools_threads: Option<usize>,
    pub needletail_parallel: bool,
    pub min_query_length: u64,
    pub min_query_coverage: f64,
    pub min_mapq: u8,
    pub preset: Option<Preset>,
    pub alignment_format: Option<AlignmentFormat>,
    pub command: Option<String>
}

/// Builder for constructing a `Scrubby` instance.
pub struct ScrubbyBuilder {
    pub input: Vec<PathBuf>,
    pub output: Vec<PathBuf>,
    pub json: Option<PathBuf>,
    pub workdir: Option<PathBuf>,
    pub read_ids: Option<PathBuf>,
    pub extract: bool,
    pub keep: bool,
    pub threads: usize,
    pub config: ScrubbyConfig,
}

impl ScrubbyBuilder {
    /// Creates a new `ScrubbyBuilder` instance with default configuration. If both aligner and classifier are `None`,
    /// the default aligner (and preset) is selected depending on file input and output length: 
    /// 
    /// * 1 - long reads with `Aligner::Minimap2` and `Preset::MapOnt`
    /// * 2 - short reads with `Aligner::Bowtie2`
    /// 
    /// If compiled with the `mm2` feature, the following default aligners are selected:
    /// 
    /// * 1 - long reads with `Aligner::Minimap2` and `Preset::MapOnt`
    /// * 2 - short reads with `Aligner::Minimap2` and `Preset::Sr`
    ///
    /// # Arguments
    ///
    /// * `input` - A vector of input file paths.
    /// * `output` - A vector of output file paths.
    ///
    /// # Example
    ///
    /// ```
    /// use scrubby::ScrubbyBuilder;
    /// use std::path::PathBuf;
    ///
    /// let builder = ScrubbyBuilder::new(
    ///     "input.fastq", 
    ///     "output.fastq"
    /// );
    /// ```
    pub fn new<P: IntoVecPathBuf>(input: P, output: P) -> Self {
        
        let input = input.into_vec_path_buf();
        let output = output.into_vec_path_buf();

        let paired_end = input.len() == 2;

        Self {
            input,
            output,
            json: None,
            workdir: None,
            read_ids: None,
            extract: false,
            keep: false,
            threads: 4,
            config: ScrubbyConfig {
                aligner: None,
                classifier: None,
                index: None,
                aligner_index: None,
                alignment: None,
                classifier_index: None,
                reads: None,
                report: None,
                taxa: Vec::new(),
                taxa_direct: Vec::new(),
                aligner_args: None,
                classifier_args: None,
                unpaired: false,
                samtools_threads: None,
                paired_end,
                needletail_parallel: true,
                min_query_length: 0,
                min_query_coverage: 0.0,
                min_mapq: 0,
                alignment_format: None,
                preset: None,
                command: None
            },
        }
    }
    /// Sets the `read_ids` field.
    ///
    /// # Example
    ///
    /// ```
    /// use scrubby::ScrubbyBuilder;
    /// use std::path::PathBuf;
    ///
    /// let builder = ScrubbyBuilder::new(...).read_ids(PathBuf::from("depleted_reads.fastq"));
    /// ```
    pub fn read_ids<T: Into<Option<PathBuf>>>(mut self, read_ids: T) -> Self {
        self.read_ids = read_ids.into();
        self
    }
    /// Sets the `json` field.
    ///
    /// # Example
    ///
    /// ```
    /// use scrubby::ScrubbyBuilder;
    /// use std::path::PathBuf;
    ///
    /// let builder = ScrubbyBuilder::new(...).json(PathBuf::from("report.json"));
    /// ```
    pub fn json<T: Into<Option<PathBuf>>>(mut self, json: T) -> Self {
        self.json = json.into();
        self
    }
    /// Sets the `command` field.
    ///
    /// # Example
    ///
    /// ```
    /// use scrubby::ScrubbyBuilder;
    /// use std::path::PathBuf;
    ///
    /// let builder = ScrubbyBuilder::new(...).command("scrubby clean ...");
    /// ```
    pub fn command<T: Into<Option<String>>>(mut self, command: T) -> Self {
        self.config.command = command.into();
        self
    }
    /// Sets the `workdir` field.
    ///
    /// # Example
    ///
    /// ```
    /// use scrubby::ScrubbyBuilder;
    /// use std::path::PathBuf;
    ///
    /// let builder = ScrubbyBuilder::new(...).workdir(PathBuf::from("workdir"));
    /// ```
    pub fn workdir<T: Into<Option<PathBuf>>>(mut self, workdir: T) -> Self {
        self.workdir = workdir.into();
        self
    }
    /// Sets the `extract` field.
    ///
    /// # Example
    ///
    /// ```
    /// use scrubby::ScrubbyBuilder;
    ///
    /// let builder = ScrubbyBuilder::new(...).extract(true);
    /// ```
    pub fn extract(mut self, extract: bool) -> Self {
        self.extract = extract;
        self
    }
    /// Sets the `keep` field.
    ///
    /// # Example
    ///
    /// ```
    /// use scrubby::ScrubbyBuilder;
    ///
    /// let builder = ScrubbyBuilder::new(...).keep(true);
    /// ```
    pub fn keep(mut self, keep: bool) -> Self {
        self.keep = keep;
        self
    }
    /// Sets the `unpaired` field.
    ///
    /// # Example
    ///
    /// ```
    /// use scrubby::ScrubbyBuilder;
    ///
    /// let builder = ScrubbyBuilder::new(...).unpaired(true);
    /// ```
    pub fn unpaired(mut self, unpaired: bool) -> Self {
        self.config.unpaired = unpaired;
        self
    }
    /// Sets the number of `threads`.
    ///
    /// # Example
    ///
    /// ```
    /// use scrubby::ScrubbyBuilder;
    ///
    /// let builder = ScrubbyBuilder::new(...).threads(8);
    /// ```
    pub fn threads(mut self, threads: usize) -> Self {
        self.threads = threads;
        self
    }
    /// Sets the `aligner` field.
    ///
    /// # Example
    ///
    /// ```
    /// use scrubby::{ScrubbyBuilder, Aligner};
    ///
    /// let builder = ScrubbyBuilder::new(...).aligner(Aligner::Minimap2);
    /// ```
    pub fn aligner<T: Into<Option<Aligner>>>(mut self, aligner: T) -> Self {
        self.config.aligner = aligner.into();
        self
    }
    /// Sets the `alignment` field.
    ///
    /// # Example
    ///
    /// ```
    /// use scrubby::ScrubbyBuilder;
    /// use std::path::PathBuf;
    ///
    /// let builder = ScrubbyBuilder::new(...).alignment(PathBuf::from("alignment.paf"));
    /// ```
    pub fn alignment<T: Into<Option<PathBuf>>>(mut self, alignment: T) -> Self {
        self.config.alignment = alignment.into();
        self
    }
    /// Sets the `alignment_format` field.
    ///
    /// # Example
    ///
    /// ```
    /// use scrubby::{ScrubbyBuilder, AlignmentFormat};
    /// use std::path::PathBuf;
    ///
    /// let builder = ScrubbyBuilder::new(...).alignment_format(AlignmentFormat::Bam);
    /// ```
    pub fn alignment_format<T: Into<Option<AlignmentFormat>>>(mut self, alignment_format: T) -> Self {
        self.config.alignment_format = alignment_format.into();
        self
    }
    /// Sets the `min_query_length` field.
    ///
    /// # Example
    ///
    /// ```
    /// use scrubby::ScrubbyBuilder;
    /// use std::path::PathBuf;
    ///
    /// let builder = ScrubbyBuilder::new(...).min_query_length(50);
    /// ```
    pub fn min_query_length(mut self, min_query_length: u64) -> Self {
        self.config.min_query_length = min_query_length;
        self
    }
    /// Sets the `min_query_coverage` field.
    ///
    /// # Example
    ///
    /// ```
    /// use scrubby::ScrubbyBuilder;
    /// use std::path::PathBuf;
    ///
    /// let builder = ScrubbyBuilder::new(...).min_query_coverage(0.5);
    /// ```
    pub fn min_query_coverage(mut self, min_query_coverage: f64) -> Self {
        self.config.min_query_coverage = min_query_coverage;
        self
    }
    /// Sets the `min_mapq` field.
    ///
    /// # Example
    ///
    /// ```
    /// use scrubby::ScrubbyBuilder;
    /// use std::path::PathBuf;
    ///
    /// let builder = ScrubbyBuilder::new(...).min_mapq(50);
    /// ```
    pub fn min_mapq(mut self, min_mapq: u8) -> Self {
        self.config.min_mapq = min_mapq;
        self
    }
    /// Sets the `classifier` field.
    ///
    /// # Example
    ///
    /// ```
    /// use scrubby::{ScrubbyBuilder, Classifier};
    ///
    /// let builder = ScrubbyBuilder::new(...).classifier(Classifier::Kraken2);
    /// ```
    pub fn classifier<T: Into<Option<Classifier>>>(mut self, classifier: T) -> Self {
        self.config.classifier = classifier.into();
        self
    }
    /// Sets the `reads` field.
    ///
    /// # Example
    ///
    /// ```
    /// use scrubby::ScrubbyBuilder;
    /// use std::path::PathBuf;
    ///
    /// let builder = ScrubbyBuilder::new(...).reads(PathBuf::from("read_classifications.tsv"));
    /// ```
    pub fn reads<T: Into<Option<PathBuf>>>(mut self, reads: T) -> Self {
        self.config.reads = reads.into();
        self
    }
    /// Sets the `report` field.
    ///
    /// # Example
    ///
    /// ```
    /// use scrubby::ScrubbyBuilder;
    /// use std::path::PathBuf;
    ///
    /// let builder = ScrubbyBuilder::new(...).report(PathBuf::from("classifier_report.tsv"));
    /// ```
    pub fn report<T: Into<Option<PathBuf>>>(mut self, report: T) -> Self {
        self.config.report = report.into();
        self
    }
    /// Sets the `index` field.
    ///
    /// # Example
    ///
    /// ```
    /// use scrubby::ScrubbyBuilder;
    /// use std::path::PathBuf;
    ///
    /// let builder = ScrubbyBuilder::new(...).index(PathBuf::from("index"));
    /// ```
    pub fn index<T: Into<Option<PathBuf>>>(mut self, index: T) -> Self {
        self.config.index = index.into();
        self
    }
    /// Sets the `aligner_index` field directly.
    ///
    /// # Example
    ///
    /// ```
    /// use scrubby::ScrubbyBuilder;
    /// use std::path::PathBuf;
    ///
    /// let builder = ScrubbyBuilder::new(...).aligner_index(PathBuf::from("aligner_index"));
    /// ```
    pub fn aligner_index<T: Into<Option<PathBuf>>>(mut self, aligner_index: T) -> Self {
        self.config.aligner_index = aligner_index.into();
        self
    }
    /// Sets the `classifier_index` field directly.
    ///
    /// # Example
    ///
    /// ```
    /// use scrubby::ScrubbyBuilder;
    /// use std::path::PathBuf;
    ///
    /// let builder = ScrubbyBuilder::new(...).classifier_index(PathBuf::from("classifier_index"));
    /// ```
    pub fn classifier_index<T: Into<Option<PathBuf>>>(mut self, classifier_index: T) -> Self {
        self.config.classifier_index = classifier_index.into();
        self
    }
    /// Sets the `taxa` field.
    ///
    /// # Example
    ///
    /// ```
    /// use scrubby::ScrubbyBuilder;
    ///
    /// let builder = ScrubbyBuilder::new(...).taxa(vec!["taxon1", "taxon2"]);
    /// ```
    pub fn taxa<T: Into<Vec<String>>>(mut self, taxa: T) -> Self {
        self.config.taxa = taxa.into();
        self
    }
    /// Sets the `taxa_direct` field.
    ///
    /// # Example
    ///
    /// ```
    /// use scrubby::ScrubbyBuilder;
    ///
    /// let builder = ScrubbyBuilder::new(...).taxa_direct(vec!["taxon_direct1", "taxon_direct2"]);
    /// ```
    pub fn taxa_direct<T: Into<Vec<String>>>(mut self, taxa_direct: T) -> Self {
        self.config.taxa_direct = taxa_direct.into();
        self
    }
    /// Sets the `classifier_args` field.
    ///
    /// # Example
    ///
    /// ```
    /// use scrubby::ScrubbyBuilder;
    ///
    /// let builder = ScrubbyBuilder::new(...).classifier_args("classifier_args".to_string());
    /// ```
    pub fn classifier_args<T: Into<Option<String>>>(mut self, classifier_args: T) -> Self {
        self.config.classifier_args = classifier_args.into();
        self
    }
    /// Sets the `aligner_args` field.
    ///
    /// # Example
    ///
    /// ```
    /// use scrubby::ScrubbyBuilder;
    ///
    /// let builder = ScrubbyBuilder::new(...).aligner_args("aligner_args".to_string());
    /// ```
    pub fn aligner_args<T: Into<Option<String>>>(mut self, aligner_args: T) -> Self {
        self.config.aligner_args = aligner_args.into();
        self
    }
    /// Sets the number of `samtools_threads`.
    ///
    /// # Example
    ///
    /// ```
    /// use scrubby::ScrubbyBuilder;
    ///
    /// let builder = ScrubbyBuilder::new(...).samtools_threads(4);
    /// ```
    pub fn samtools_threads<T: Into<Option<usize>>>(mut self, threads: T) -> Self {
        self.config.samtools_threads = threads.into();
        self
    }
    /// Sets the `needletail_parallel` field.
    ///
    /// # Example
    ///
    /// ```
    /// use scrubby::ScrubbyBuilder;
    ///
    /// let builder = ScrubbyBuilder::new(...).needletail_parallel(true);
    /// ```
    pub fn needletail_parallel(mut self, parallel: bool) -> Self {
        self.config.needletail_parallel = parallel;
        self
    }
    /// Sets the minimap2 preset `preset` field.
    ///
    /// # Example
    ///
    /// ```
    /// use scrubby::ScrubbyBuilder;
    /// use scrubby::Preset;
    ///
    /// let builder = ScrubbyBuilder::new(...).preset(Preset::MapOnt);
    /// ```
    pub fn preset<T: Into<Option<Preset>>>(mut self, preset: T) -> Self {
        self.config.preset = preset.into();
        self
    }
    pub fn validate_base_config(&mut self) -> Result<(), ScrubbyError> {

        // Check if input and output vectors are not empty
        if self.input.is_empty() || self.output.is_empty() {
            return Err(ScrubbyError::EmptyInputOutput);
        }
        // Check if input and output vectors have the same length
        if self.input.len() != self.output.len() {
            return Err(ScrubbyError::MismatchedInputOutputLength);
        }
        // Check if input and output vectors length is limited to one or two
        if self.input.len() > 2 || self.output.len() > 2 {
            return Err(ScrubbyError::InputOutputLengthExceeded);
        }
        // Check if each input file exists and is a file
        for input_file in &self.input {
            if !input_file.exists() || !input_file.is_file() {
                return Err(ScrubbyError::MissingInputReadFile(input_file.clone()));
            }
        }
        // If a workdir is provided, check if it exists and is a directory, otherwise create it
        if let Some(dir) = &self.workdir {
            if !dir.exists() || !dir.is_dir() {
                create_dir_all(&dir)?;
            }
        }

        if self.config.index.is_some() {
            if self.config.aligner.is_some() {
                self.config.aligner_index = self.config.index.clone()
            } else if self.config.classifier.is_some() {
                self.config.classifier_index = self.config.index.clone()
            } else {
                // Default aligner index and settings
                self.config.aligner_index = self.config.index.clone()
            }
        }

        Ok(())
    }
    /// Builds the `Scrubby` instance.
    ///
    /// # Returns
    ///
    /// * `Result<Scrubby, ScrubbyError>` - Ok with the constructed Scrubby instance, otherwise an error.
    ///
    /// # Example
    ///
    /// ```
    /// use scrubby::ScrubbyBuilder;
    /// 
    /// let scrubby = ScrubbyBuilder::new(...).build().unwrap();
    /// ```
    pub fn build(mut self) -> Result<Scrubby, ScrubbyError> {

        self.validate_base_config()?;

        // Check if either aligner or classifier is set
        if self.config.aligner.is_none() && self.config.classifier.is_none() {

            // Default aligners depend on feature configuration
            #[cfg(not(feature = "mm2"))]
            {
                self.config.aligner = if self.config.paired_end { Some(Aligner::Bowtie2) } else { Some(Aligner::Minimap2) }
            }
            #[cfg(feature = "mm2")]
            {
                self.config.aligner = Some(Aligner::Minimap2Rs)
            }
        }
        // Check if only one of aligner or classifier is set
        if self.config.aligner.is_some() && self.config.classifier.is_some() {
            return Err(ScrubbyError::AlignerAndClassifierConfigured);
        }
        // Check if only one of aligner or classifier index is set
        if self.config.aligner_index.is_some() && self.config.classifier_index.is_some() {
            return Err(ScrubbyError::AlignerAndClassifierIndexConfigured);
        }
        // Check if classifier is set and necessary fields are populated
        if let Some(_) = &self.config.classifier {
            if self.config.classifier_index.is_none() {
                return Err(ScrubbyError::MissingClassifierIndex);
            }
            if self.config.taxa.is_empty() && self.config.taxa_direct.is_empty() {
                return Err(ScrubbyError::MissingTaxa);
            }
        }
        // Check if aligner is set and necessary fields are populated
        if let Some(_) = &self.config.aligner {
            if self.config.aligner_index.is_none() {
                return Err(ScrubbyError::MissingAlignmentIndex);
            }
        }
        // Check if the classifier index exists and is a directory
        if let Some(dir) = &self.config.classifier_index {
            if !dir.exists() || !dir.is_dir() {
                return Err(ScrubbyError::MissingClassifierIndexDirectory(dir.clone()));
            }
        }
        
        // If the index file for Strobealign ends in ".sti" strobealign expects the 
        // underlying FASTA file to be in the same directory (v0.13.0) - this is 
        // kinda weird...
        if let Some(Aligner::Strobealign) = &self.config.aligner {
            if let Some(file) = &self.config.aligner_index {
                if file.extension().unwrap_or_default() == "sti" {
                    let index_base_file = file.with_extension("").with_extension("");
                    if !index_base_file.exists() {
                        return Err(ScrubbyError::MissingStrobealignIndexBaseFile(index_base_file.clone()));
                    }
                }
            }
        }
        // If Bowtie2 aligner is set, check index files exist and are files
        // otherwise check if the aligner index file provided exists and is a file
        if let Some(Aligner::Bowtie2) = &self.config.aligner {
            if let Some(file) = &self.config.aligner_index {
                // Check if Bowtie2 index files are all present
                let bowtie2_small_extensions = ["1.bt2", "2.bt2", "3.bt2", "4.bt2", "rev.1.bt2", "rev.2.bt2"];
                let bowtie2_large_extensions = ["1.bt21", "2.bt21", "3.bt21", "4.bt21", "rev.1.bt21", "rev.2.bt21"];
                for (small_ext, large_ext) in bowtie2_small_extensions.iter().zip(bowtie2_large_extensions.iter()) {
                    let small_index_file = file.with_extension(small_ext);
                    let large_index_file = file.with_extension(large_ext);
                    if !small_index_file.exists() || !small_index_file.is_file() {
                        if !large_index_file.exists() || !large_index_file.is_file() {
                            return Err(ScrubbyError::MissingBowtie2IndexFiles(file.clone()));
                        }
                    }
                }
            }
        } else {
            if let Some(file) = &self.config.aligner_index {
                if !file.exists() || !file.is_file() {
                    return Err(ScrubbyError::MissingAlignmentIndexFile(file.clone()));
                }
            }
        }



        // Check that a default preset is set with Minimap2 or that 
        // the preset is supported by Minimap2
        if let Some(Aligner::Minimap2) = &self.config.aligner {
            match self.config.preset {
                None => {
                    if self.config.paired_end {
                        self.config.preset = Some(Preset::Sr)
                    } else {
                        self.config.preset = Some(Preset::MapOnt)
                    }
                },
                Some(ref preset) => {
                    if [Preset::Lr].contains(preset) {
                        return Err(ScrubbyError::Minimap2PresetNotSupported(preset.to_owned()))
                    }
                }
            }
        }
        // Check that a default preset is set with Minigraph
        if let Some(Aligner::Minigraph) = &self.config.aligner {
            match self.config.preset {
                None => {
                    if self.config.paired_end {
                        self.config.preset = Some(Preset::Sr)
                    } else {
                        self.config.preset = Some(Preset::Lr)
                    }
                },
                Some(ref preset) => {
                    if ![Preset::Lr, Preset::Sr, Preset::Asm].contains(preset) {
                        return Err(ScrubbyError::MinigraphPresetNotSupported(preset.to_owned()))
                    }
                }
            }
        }
        #[cfg(feature = "mm2")]
        if let Some(Aligner::Minimap2Rs) = &self.config.aligner {
            match self.config.preset {
                None => {
                    if self.config.paired_end {
                        self.config.preset = Some(Preset::Sr)
                    } else {
                        self.config.preset = Some(Preset::MapOnt)
                    }
                },
                Some(ref preset) => {
                    if [Preset::Lr].contains(preset) {
                        return Err(ScrubbyError::Minimap2PresetNotSupported(preset.to_owned()))
                    }
                }
            }
        }

        Ok(Scrubby {
            input: self.input,
            output: self.output,
            read_ids: self.read_ids,
            json: self.json,
            workdir: self.workdir,
            extract: self.extract,
            keep: self.keep,
            threads: self.threads,
            config: self.config,
        })
    }
    /// Builds the `Scrubby` instance with the classifier output cleaning configuration.
    ///
    /// # Returns
    ///
    /// * `Result<Scrubby, ScrubbyError>` - Ok with the constructed Scrubby instance, otherwise an error.
    ///
    /// # Example
    ///
    /// ```
    /// use scrubby::ScrubbyBuilder;
    /// 
    /// let scrubby = ScrubbyBuilder::new(...).build_classifier().unwrap();
    /// ```
    pub fn build_classifier(mut self) -> Result<Scrubby, ScrubbyError> {

        self.validate_base_config()?;

        // Check if classifier read classifications is set
        if self.config.reads.is_none() {
            return Err(ScrubbyError::MissingClassifierReadClassfications);
        }
        // Check if classifier read classifications report is set
        if self.config.report.is_none() {
            return Err(ScrubbyError::MissingClassifierClassificationReport);
        }
        // Check if taxa directive for cleaning is set
        if self.config.taxa.is_empty() && self.config.taxa_direct.is_empty() {
            return Err(ScrubbyError::MissingTaxa);
        }

        Ok(Scrubby {
            input: self.input,
            output: self.output,
            read_ids: self.read_ids,
            json: self.json,
            workdir: self.workdir,
            extract: self.extract,
            keep: self.keep,
            threads: self.threads,
            config: self.config,
        })
    }/// Builds the `Scrubby` instance with the alignment output cleaning configuration.
    ///
    /// # Returns
    ///
    /// * `Result<Scrubby, ScrubbyError>` - Ok with the constructed Scrubby instance, otherwise an error.
    ///
    /// # Example
    ///
    /// ```
    /// use scrubby::ScrubbyBuilder;
    /// 
    /// let scrubby = ScrubbyBuilder::new(...).build_alignment().unwrap();
    /// ```
    pub fn build_alignment(mut self) -> Result<Scrubby, ScrubbyError> {

        self.validate_base_config()?;

        if self.config.alignment.is_none() {
            return Err(ScrubbyError::MissingAlignment);
        }

        Ok(Scrubby {
            input: self.input,
            output: self.output,
            read_ids: self.read_ids,
            json: self.json,
            workdir: self.workdir,
            extract: self.extract,
            keep: self.keep,
            threads: self.threads,
            config: self.config,
        })
    }
}
