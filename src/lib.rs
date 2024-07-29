pub mod scrubby;
pub mod error;
pub mod utils;
pub mod terminal;
pub mod cleaner;
pub mod classifier;
pub mod alignment;
pub mod download;
pub mod report;

pub mod prelude {
    pub use crate::scrubby::{Aligner, Classifier, Scrubby, ScrubbyBuilder};
    pub use crate::download::{ScrubbyDownloader, ScrubbyDownloaderBuilder, ScrubbyIndex};
    pub use crate::alignment::AlignmentFormat;
    pub use crate::error::ScrubbyError;
    pub use crate::utils::init_logger;
}
