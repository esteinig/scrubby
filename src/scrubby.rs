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
//! The primary structures and enumerations provided are `Aligner`, `Classifier`, 
//! `Scrubby`, `ScrubbyConfig`, and `ScrubbyBuilder`. These components can be imported
//! through the prelude module: `scrubby::prelude::*`.

use serde::{Serialize, Deserialize};
use std::fs::create_dir_all;
use std::path::PathBuf;
use std::fmt;
use crate::cleaner::Cleaner;
use crate::error::ScrubbyError;

/// Enum representing the available aligners.
#[derive(Serialize, Deserialize, Clone, Debug, clap::ValueEnum)]
pub enum Aligner {
    Bowtie2,
    Minimap2,
    Strobealign,
}

impl fmt::Display for Aligner {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Aligner::Bowtie2 => write!(f, "Bowtie2"),
            Aligner::Minimap2 => write!(f, "Minimap2"),
            Aligner::Strobealign => write!(f, "Strobealign"),
        }
    }
}

/// Enum representing the available classifiers.
#[derive(Serialize, Deserialize, Clone, Debug, clap::ValueEnum)]
pub enum Classifier {
    Kraken2,
    Metabuli,
}

impl fmt::Display for Classifier {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Classifier::Kraken2 => write!(f, "Kraken2"),
            Classifier::Metabuli => write!(f, "Metabuli"),
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
    pub reverse: bool,
    pub keep: bool,
    pub threads: usize,
    pub config: ScrubbyConfig,
}

impl Scrubby {
    /// Creates a new `Scrubby` instance using the provided parameters.
    ///
    /// # Arguments
    ///
    /// * `input` - A vector of input file paths.
    /// * `output` - A vector of output file paths.
    /// * `aligner` - Optional aligner configuration.
    /// * `aligner_index` - Optional path to the aligner index.
    /// * `classifier` - Optional classifier configuration.
    /// * `classifier_index` - Optional path to the classifier index.
    ///
    /// # Example
    ///
    /// ```
    /// use scrubby::{Scrubby, Aligner, Classifier};
    /// use std::path::PathBuf;
    ///
    /// let scrubby = Scrubby::new(
    ///     vec![PathBuf::from("input.fastq")], 
    ///      vec![PathBuf::from("output.fastq")], 
    ///     Some(Aligner::Bowtie2), 
    ///     None, 
    ///     Some(Classifier::Kraken2), 
    ///     Some(PathBuf::from("kraken2_db/")
    /// )).unwrap();
    /// ```
    pub fn new(
        input: Vec<PathBuf>, 
        output: Vec<PathBuf>, 
        aligner: Option<Aligner>, 
        aligner_index: Option<PathBuf>, 
        classifier: Option<Classifier>, 
        classifier_index: Option<PathBuf>
    ) -> Result<Self, ScrubbyError> {
        ScrubbyBuilder::new(input, output)
            .aligner(aligner)
            .classifier(classifier)
            .aligner_index(aligner_index)
            .classifier_index(classifier_index)
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
    ///     vec![PathBuf::from("input.fastq")], 
    ///     vec![PathBuf::from("output.fastq")]
    /// );
    /// ```
    pub fn builder(input: Vec<PathBuf>, output: Vec<PathBuf>) -> ScrubbyBuilder {
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

        Ok(())
    }
}

/// Configuration structure for the Scrubby tool.
#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct ScrubbyConfig {
    pub aligner: Option<Aligner>,
    pub classifier: Option<Classifier>,
    pub aligner_index: Option<PathBuf>,
    pub classifier_index: Option<PathBuf>,
    pub taxa: Vec<String>,
    pub taxa_direct: Vec<String>,
    pub classifier_args: Option<String>,
    pub aligner_args: Option<String>,
    pub unpaired: bool,
    pub paired_end: bool,
    pub samtools_threads: Option<usize>,
    pub needletail_parallel: bool,
}

/// Builder for constructing a `Scrubby` instance.
pub struct ScrubbyBuilder {
    pub input: Vec<PathBuf>,
    pub output: Vec<PathBuf>,
    pub json: Option<PathBuf>,
    pub workdir: Option<PathBuf>,
    pub read_ids: Option<PathBuf>,
    pub reverse: bool,
    pub keep: bool,
    pub threads: usize,
    pub config: ScrubbyConfig,
}

impl ScrubbyBuilder {
    /// Creates a new `ScrubbyBuilder` instance.
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
    ///     vec![PathBuf::from("input.fastq")], 
    ///     vec![PathBuf::from("output.fastq")]
    /// );
    /// ```
    pub fn new(input: Vec<PathBuf>, output: Vec<PathBuf>) -> Self {
        
        let paired_end = input.len() == 2;

        Self {
            input,
            output,
            json: None,
            workdir: None,
            read_ids: None,
            reverse: false,
            keep: false,
            threads: 4,
            config: ScrubbyConfig {
                aligner: None,
                classifier: None,
                aligner_index: None,
                classifier_index: None,
                taxa: Vec::new(),
                taxa_direct: Vec::new(),
                aligner_args: None,
                classifier_args: None,
                unpaired: false,
                samtools_threads: None,
                paired_end,
                needletail_parallel: true,
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
    /// let builder = ScrubbyBuilder::new(...).read_ids(PathBuf::from("reads.fastq"));
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
    /// let builder = ScrubbyBuilder::new(...).json(PathBuf::from("config.json"));
    /// ```
    pub fn json<T: Into<Option<PathBuf>>>(mut self, json: T) -> Self {
        self.json = json.into();
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
    /// Sets the `reverse` field.
    ///
    /// # Example
    ///
    /// ```
    /// use scrubby::ScrubbyBuilder;
    ///
    /// let builder = ScrubbyBuilder::new(...).reverse(true);
    /// ```
    pub fn reverse(mut self, reverse: bool) -> Self {
        self.reverse = reverse;
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
    /// Sets the `aligner_index` field.
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
    /// Sets the `classifier_index` field.
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
    /// let builder = ScrubbyBuilder::new(...).taxa(vec!["taxon1".to_string(), "taxon2".to_string()]);
    /// ```
    pub fn taxa(mut self, taxa: Vec<String>) -> Self {
        self.config.taxa = taxa;
        self
    }
    /// Sets the `taxa_direct` field.
    ///
    /// # Example
    ///
    /// ```
    /// use scrubby::ScrubbyBuilder;
    ///
    /// let builder = ScrubbyBuilder::new(...).taxa_direct(vec!["taxon_direct1".to_string(), "taxon_direct2".to_string()]);
    /// ```
    pub fn taxa_direct(mut self, taxa_direct: Vec<String>) -> Self {
        self.config.taxa_direct = taxa_direct;
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
    pub fn build(self) -> Result<Scrubby, ScrubbyError> {

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
        // Check if either aligner or classifier is set
        if self.config.aligner.is_none() && self.config.classifier.is_none() {
            return Err(ScrubbyError::MissingClassifierOrAligner);
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

        Ok(Scrubby {
            input: self.input,
            output: self.output,
            read_ids: self.read_ids,
            json: self.json,
            workdir: self.workdir,
            reverse: self.reverse,
            keep: self.keep,
            threads: self.threads,
            config: self.config,
        })
    }
}
