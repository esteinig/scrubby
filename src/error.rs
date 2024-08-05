use thiserror::Error;
use std::path::PathBuf;
use crate::scrubby::{Aligner, Classifier, Preset};

/// Represents all possible errors that can occur in the Scrubby application.
#[derive(Error, Debug)]
pub enum ScrubbyError {
    #[cfg(feature = "htslib")]
    /// Indicates failure to parse a BAM file
    #[error("failed to parse records from BAM")]
    HtslibError(#[from] rust_htslib::errors::Error),
    /// Represents all other cases of `std::io::Error`.
    #[error(transparent)]
    IoError(#[from] std::io::Error),
    /// Represents errors from building a Rayon thread pool.
    #[error(transparent)]
    RayonThreadPoolError(#[from] rayon::ThreadPoolBuildError),
    /// Represents all other cases of `niffler::Error`.
    #[error(transparent)]
    NifflerError(#[from] niffler::Error),
    /// Represents all other cases of `needletail::errors::ParseError`.
    #[error(transparent)]
    NeedletailParseError(#[from] needletail::errors::ParseError),
    /// Represents all other cases of `reqwest::Error`.
    #[error(transparent)]
    ReqwestError(#[from] reqwest::Error),
    /// Represents all other cases of `csv::Error`.
    #[error(transparent)]
    CsvError(#[from] csv::Error),
    /// Represents all other cases of `serde_json::Error`.
    #[error(transparent)]
    SerdeJsonError(#[from] serde_json::Error),
    /// Failed to make the download request
    #[error("failed to execute request: {0}")]
    DownloadFailedRequest(reqwest::StatusCode),
    /// Failed to configure the downloader through the builder pattern due to missing field
    #[error("failed to configure the output directory field for the downloader")]
    DownloaderMissingOutdir,
    /// Indicates failure to parse a record name from BAM file
    #[error("failed to parse record name from BAM")]
    RecordNameUtf8Error(#[from] std::str::Utf8Error),
    /// Indicates failure to parse a target name from BAM file
    #[error("failed to parse a valid record target name from BAM")]
    RecordTargetIdError(#[from] std::num::TryFromIntError),
    /// Indicates failure to parse an u64 from PAF
    #[error("failed to parse a valid integer from PAF")]
    PafRecordIntegerError(#[from] std::num::ParseIntError),
    /// Represents an error when failing to extract a sequence record header.
    #[error("failed to extract sequence record header")]
    NeedletailHeader(#[source] std::str::Utf8Error),
    /// Represents an error when failing to extract a valid header of a read.
    #[error("failed to extract a valid header of read")]
    NeedletailFastqHeader,
    /// Represents an error when both aligner and classifier are configured simultaneously.
    #[error("Unable to specify both aligner and classifier.")]
    AlignerAndClassifierConfigured,
    /// Represents an error when both aligner and classifier indices are specified simultaneously.
    #[error("Unable to specify both aligner and classifier indices.")]
    AlignerAndClassifierIndexConfigured,
    /// Represents an error when the alignment format is not explicitly set and not recognized from extension
    #[error("Unable to recognize alignment input format from extension.")]
    AlignmentInputFormatNotRecognized,
    /// Represents an error when the alignment format is explicitly set and not recognized
    #[error("Unable to recognize alignment input format - is this version compiled with 'htslib'?")]
    AlignmentInputFormatInvalid,
    /// Represents an error when input and output lengths do not match.
    #[error("Input and output must be of the same length.")]
    MismatchedInputOutputLength,
    /// Represents an error when classifier is set but `taxa` or `taxa_direct` is empty.
    #[error("If classifier is set, `taxa` or `taxa_direct` must not be empty.")]
    MissingTaxa,
    /// Represents an error when classifier index is not set while classifier is configured.
    #[error("Classifier index must be set when classifier is configured.")]
    MissingClassifierIndex,
    /// Represents an error when classifier read classfication file is not set while classifier cleaning procedure is configured.
    #[error("Classifier read classification input must be set when classifier cleaning procedure is configured.")]
    MissingClassifierReadClassfications,
    /// Represents an error when classifier read classfication report is not set while classifier cleaning procedure is configured.
    #[error("Classifier read classification report input must be set when classifier cleaning procedure is configured.")]
    MissingClassifierClassificationReport,
    /// Represents an error when alignment index is not set while aligner is configured.
    #[error("Alignment index must be set when aligner is configured.")]
    MissingAlignmentIndex,
    /// Represents an error when alignment output is not set while alignment is configured.
    #[error("Alignment output must be set when alignment is configured.")]
    MissingAlignment,
    /// Represents an error when neither classifier nor aligner is set.
    #[error("Either classifier or aligner must be set.")]
    MissingClassifierOrAligner,
    /// Represents an error when input and output vectors are empty.
    #[error("Input and output vectors must not be empty.")]
    EmptyInputOutput,
    /// Represents an error when input and output vectors contain more than two elements.
    #[error("Input and output vectors must not contain more than two elements.")]
    InputOutputLengthExceeded,
    /// Represents an error when a command execution fails.
    #[error("Failed to execute command '{0}': {1}")]
    CommandExecutionFailed(String, String),
    /// Represents an error when a command exits with a non-zero status code.
    #[error("Command '{0}' exited with status code: {1}")]
    CommandFailed(String, i32),
    /// Represents an error when no aligner is configured.
    #[error("No aligner configured.")]
    MissingAligner,
    /// Represents an error when no classifier is configured.
    #[error("No classifier configured.")]
    MissingClassifier,
    /// Represents an error when no preset is configured.
    #[error("Minimap2 was set as aligner but no preset was configured.")]
    MissingMinimap2Preset,
    /// Represents an error when no preset is configured.
    #[error("Minigraph was set as aligner but no preset was configured.")]
    MissingMinigraphPreset,
    /// Represents an error when the strobealign index base file is not found.
    #[error("Strobealign index file provided but matching base file was not found in the same directory (required): {0}")]
    MissingStrobealignIndexBaseFile(PathBuf),
    /// Represents an error when the input read file is not found.
    #[error("Read input file was not found: {0}")]
    MissingInputReadFile(PathBuf),
    /// Represents an error when the alignment index file is not found.
    #[error("Alignment index file was not found: {0}")]
    MissingAlignmentIndexFile(PathBuf),
    /// Represents an error when neither small nor large index files for Bowtie2 are found with the specified base path.
    #[error("Neither small nor large index files for Bowtie2 were found with base path: {0}")]
    MissingBowtie2IndexFiles(PathBuf),
    /// Represents an error when the classifier index directory is not found.
    #[error("Classifier index directory was not found: {0}")]
    MissingClassifierIndexDirectory(PathBuf),
    /// Represents an error when the specified aligner cannot be executed, possibly due to it not being installed.
    #[error("Aligner `{0}` cannot be executed - is it installed?")]
    AlignerDependencyMissing(Aligner),
    /// Represents an error when the specified classifier cannot be executed, possibly due to it not being installed.
    #[error("Classifier `{0}` cannot be executed - is it installed?")]
    ClassifierDependencyMissing(Classifier),
    /// Represents a failure to count a taxonomic parent during report parsing from `Kraken2`.
    #[error("failed to provide a parent taxon while parsing report from `Kraken2`")]
    KrakenReportTaxonParent,
    /// Represents a failure to convert the read field from string to numeric field in the report file from `Kraken2`.
    #[error("failed to convert the read field in the report from `Kraken2`")]
    KrakenReportReadFieldConversion,
    /// Represents a failure to convert the direct read field from string to numeric field in the report file from `Kraken2`.
    #[error("failed to convert the direct read field in the report from `Kraken2`")]
    KrakenReportDirectReadFieldConversion,
    /// Represents an error when the aligner builder fails for `minimap2-rs`
    #[error("Failed to build aligner with `minimap2-rs`: {0}")]
    Minimap2RustAlignerBuilderFailed(String),
    /// Represents an error when the aligner builder fails for `minimap2-rs`
    #[error("Failed to align read with `minimap2-rs`: {0}")]
    Minimap2RustAlignmentFailed(String),
    /// Represents an error when an unsupported preset is set for `minimap2`
    #[error("Preset not supported for `minimap2` or `minimap2-rs`: {0}")]
    Minimap2PresetNotSupported(Preset),
    /// Represents an error when an unsupported preset is set for `minigraph`
    #[error("Preset not supported for `minigraph`: {0}")]
    MinigraphPresetNotSupported(Preset),
    /// Represents an error when a model save operation fails
    #[error("failed to save neural network model")]
    SaveNeuralNetworkModel,
    /// Represents an error when a model read operation fails
    #[error("failed to read neural network model")]
    ReadNeuralNetworkModel,
    /// Represents an error in label extraction function
    #[error("failed to read label from training data file; this should be a numeric suffix to the filename without extensions")]
    ReadNeuralNetworkModelLabel,
    /// Represents an error when a model read operation fails
    #[error("failed to read input sequence file: {0}")]
    ReadNeuralNetworkFastq(PathBuf),
}
