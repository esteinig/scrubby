use serde::{Serialize, Deserialize};
use std::path::PathBuf;
use std::fmt;
use crate::cleaner::Cleaner;
use crate::error::ScrubbyError;

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

#[derive(Serialize, Deserialize, Clone, Debug, clap::ValueEnum)]
pub enum Classifier {
    Kraken2,
    Metabuli,
}
impl fmt::Display for Classifier {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Classifier::Kraken2 => write!(f, "Kraken2"),
            Classifier::Metabuli => write!(f, "Metabuli")
        }
    }
}

#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct Scrubby {
    pub input: Vec<PathBuf>,
    pub output: Vec<PathBuf>,
    pub reads: Option<PathBuf>,
    pub json: Option<PathBuf>,
    pub workdir: Option<PathBuf>,
    pub reverse: bool,
    pub keep: bool,
    pub threads: usize,
    pub config: ScrubbyConfig,
}
impl Scrubby {
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
    pub fn builder(
        input: Vec<PathBuf>, 
        output: Vec<PathBuf>
    ) -> ScrubbyBuilder {
        ScrubbyBuilder::new(input, output)
    }
    pub fn clean(self) -> Result<(), ScrubbyError> {

        let cleaner = Cleaner::from_scrubby(self)?;

        cleaner.run_alignment()?;

        Ok(())
    }
}

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
    pub unpaired: bool
}

pub struct ScrubbyBuilder {
    pub input: Vec<PathBuf>,
    pub output: Vec<PathBuf>,
    pub reads: Option<PathBuf>,
    pub json: Option<PathBuf>,
    pub workdir: Option<PathBuf>,
    pub reverse: bool,
    pub keep: bool,
    pub threads: usize,
    pub config: ScrubbyConfig,
}
impl ScrubbyBuilder {
    pub fn new(input: Vec<PathBuf>, output: Vec<PathBuf>) -> Self {
        Self {
            input,
            output,
            reads: None,
            json: None,
            workdir: None,
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
                unpaired: false
            },
        }
    }
    pub fn reads<T: Into<Option<PathBuf>>>(mut self, reads: T) -> Self {
        self.reads = reads.into();
        self
    }
    pub fn json<T: Into<Option<PathBuf>>>(mut self, json: T) -> Self {
        self.json = json.into();
        self
    }
    pub fn workdir<T: Into<Option<PathBuf>>>(mut self, workdir: T) -> Self {
        self.workdir = workdir.into();
        self
    }
    pub fn reverse(mut self, reverse: bool) -> Self {
        self.reverse = reverse;
        self
    }
    pub fn keep(mut self, keep: bool) -> Self {
        self.keep = keep;
        self
    }
    pub fn unpaired(mut self, unpaired: bool) -> Self {
        self.config.unpaired = unpaired;
        self
    }
    pub fn threads(mut self, threads: usize) -> Self {
        self.threads = threads;
        self
    }
    pub fn aligner<T: Into<Option<Aligner>>>(mut self, aligner: T) -> Self {
        self.config.aligner = aligner.into();
        self
    }
    pub fn classifier<T: Into<Option<Classifier>>>(mut self, classifier: T) -> Self {
        self.config.classifier = classifier.into();
        self
    }
    pub fn aligner_index<T: Into<Option<PathBuf>>>(mut self, aligner_index: T) -> Self {
        self.config.aligner_index = aligner_index.into();
        self
    }
    pub fn classifier_index<T: Into<Option<PathBuf>>>(mut self, classifier_index: T) -> Self {
        self.config.classifier_index = classifier_index.into();
        self
    }
    pub fn taxa(mut self, taxa: Vec<String>) -> Self {
        self.config.taxa = taxa;
        self
    }
    pub fn taxa_direct(mut self, taxa_direct: Vec<String>) -> Self {
        self.config.taxa_direct = taxa_direct;
        self
    }
    pub fn classifier_args<T: Into<Option<String>>>(mut self, classifier_args: T) -> Self {
        self.config.classifier_args = classifier_args.into();
        self
    }
    pub fn aligner_args<T: Into<Option<String>>>(mut self, aligner_args: T) -> Self {
        self.config.aligner_args = aligner_args.into();
        self
    }
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

        // Check if either aligner or classifier is set - otherwise use default aligner
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
            if self.config.taxa.is_empty() || self.config.taxa_direct.is_empty() {
                return Err(ScrubbyError::MissingTaxa);
            }
        }

        // Check if aligner is set and necessary fields are populated
        if let Some(_) = &self.config.aligner {
            if self.config.aligner_index.is_none() {
                return Err(ScrubbyError::MissingAlignmentIndex);
            }
        }

        // Check if each input file exists and is a file
        for input_file in &self.input {
            if !input_file.exists() || !input_file.is_file() {
                return Err(ScrubbyError::MissingInputReadFile(input_file.clone()));
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

        if let Some(Aligner::Bowtie2) = &self.config.aligner {
            // Check if Bowtie2 index files are all present
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
            reads: self.reads,
            json: self.json,
            workdir: self.workdir,
            reverse: self.reverse,
            keep: self.keep,
            threads: self.threads,
            config: self.config,
        })
    }
}