pub mod scrubby;
pub mod error;
pub mod utils;
pub mod terminal;
pub mod cleaner;
pub mod classifier;

pub mod prelude {
    pub use crate::scrubby::{Aligner, Classifier, Scrubby};
    pub use crate::cleaner::Cleaner;
    pub use crate::error::ScrubbyError;
    pub use crate::utils::init_logger;
}
