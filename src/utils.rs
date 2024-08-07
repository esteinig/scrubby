use std::ffi::OsStr;
use std::path::{Path, PathBuf};
use crate::scrub::ScrubberError;

pub trait CompressionExt {
    fn from_path<S: AsRef<OsStr> + ?Sized>(p: &S) -> Self;
}

pub trait FileNameString {
    fn file_name_string<S: AsRef<OsStr> + ?Sized>(p: &S) -> Result<String, ScrubberError>;
}

impl FileNameString for PathBuf {
    fn file_name_string<S: AsRef<OsStr> + ?Sized>(p: &S) -> Result<String, ScrubberError> {
        let path = PathBuf::from(p);
        if let Some(file_name) = path.file_name() {
            if let Some(file_name_str) = file_name.to_str() {
                Ok(file_name_str.to_string())
            } else {
                Err(ScrubberError::FileNameInvalidUtf8)
            }
        } else {
            Err(ScrubberError::FileNameNotFound)
        }
    }
}

/// Attempts to infer the compression type from the file extension.
/// If the extension is not known, then Uncompressed is returned.
impl CompressionExt for niffler::compression::Format {
    fn from_path<S: AsRef<OsStr> + ?Sized>(p: &S) -> Self {
        let path = Path::new(p);
        match path.extension().map(|s| s.to_str()) {
            Some(Some("gz")) => Self::Gzip,
            Some(Some("bz") | Some("bz2")) => Self::Bzip,
            Some(Some("lzma")) => Self::Lzma,
            _ => Self::No,
        }
    }
}

pub fn get_file_strings_from_input(
    input: &Vec<PathBuf>,
) -> Result<[Option<String>; 2], ScrubberError> {
    match input.len() {
        2 => {
            let file1 = input[0]
                .clone()
                .into_os_string()
                .into_string()
                .map_err(|_| ScrubberError::InvalidFilePathConversion)?;
            let file2 = input[1]
                .clone()
                .into_os_string()
                .into_string()
                .map_err(|_| ScrubberError::InvalidFilePathConversion)?;
            Ok([Some(file1), Some(file2)])
        }
        1 => Ok([
            Some(
                input[0]
                    .clone()
                    .into_os_string()
                    .into_string()
                    .map_err(|_| ScrubberError::InvalidFilePathConversion)?,
            ),
            None,
        ]),
        _ => Err(ScrubberError::FileNumberError),
    }
}

