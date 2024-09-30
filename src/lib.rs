pub mod scrubby;
pub mod error;
pub mod utils;
pub mod terminal;
pub mod cleaner;
pub mod classifier;
pub mod alignment;
pub mod download;
pub mod report;

#[cfg(feature = "nn")]
pub mod identity;

pub mod prelude {
    pub use crate::download::{ScrubbyDownloader, ScrubbyDownloaderBuilder, ScrubbyIndex};
    pub use crate::scrubby::{Aligner, Classifier, Preset, Scrubby, ScrubbyConfig, ScrubbyBuilder};
    pub use crate::utils::{ReadDifference, ReadDifferenceBuilder};
    pub use crate::alignment::{ReadAlignment, AlignmentFormat};
    pub use crate::report::ScrubbyReport;
    pub use crate::error::ScrubbyError;
}
