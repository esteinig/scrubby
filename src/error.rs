use std::path::PathBuf;

use thiserror::Error;

use crate::scrubby::{Aligner, Classifier};

#[derive(Error, Debug)]
pub enum ScrubbyError {
    
    #[error("Unable to specify both aligner and classifier.")]
    AlignerAndClassifierConfigured,

    #[error("Unable to specify both aligner and classifier indices.")]
    AlignerAndClassifierIndexConfigured,

    #[error("Input and output must be of the same length.")]
    MismatchedInputOutputLength,
    
    #[error("If classifier is set, `taxa` or `taxa_direct` must not be empty.")]
    MissingTaxa,
    
    #[error("Classifier index must be set when classifier is configured.")]
    MissingClassifierIndex,
    
    #[error("Alignment index must be set when aligner is configured.")]
    MissingAlignmentIndex,
    
    #[error("Either classifier or aligner must be set.")]
    MissingClassifierOrAligner,

    #[error("Input and output vectors must not be empty.")]
    EmptyInputOutput,

    #[error("Input and output vectors must not contain more than two elements.")]
    InputOutputLengthExceeded,

    #[error("Failed to execute command '{0}': {1}")]
    CommandExecutionFailed(String, String),

    #[error("Command '{0}' exited with status code: {1}")]
    CommandFailed(String, i32),

    #[error("No aligner configured.")]
    MissingAligner,

    #[error("Strobealign index file provided but matching base file was not found in the same directory (required): {0}")]
    MissingStrobealignIndexBaseFile(PathBuf),

    #[error("Read input file was not found: {0}")]
    MissingInputReadFile(PathBuf),

    #[error("Alignment index file was not found: {0}")]
    MissingAlignmentIndexFile(PathBuf),

    #[error("Neither small nor large index files for Bowtie2 were found with base path: {0}")]
    MissingBowtie2IndexFiles(PathBuf),

    #[error("Classifier index directory was not found: {0}")]
    MissingClassifierIndexDirectory(PathBuf),

    #[error("Aligner `{0}` cannot be executed - is it installed?")]
    AlignerDependencyMissing(Aligner),

    #[error("Classifier `{0}` cannot be executed - is it installed?")]
    ClassifierDependencyMissing(Classifier),
}