use anyhow::{Result, Context};
use std::path::PathBuf;
use needletail::{parse_fastx_file, Sequence};
use std::str::from_utf8;
use std::process::Command;
use thiserror::Error;
use chrono::Local;
use std::fs;
use std::fs::{create_dir_all, canonicalize};

// See: https://nick.groenen.me/posts/rust-error-handling/

/// WordCountError enumerates all possible errors returned by this library.
#[derive(Error, Debug)]
pub enum ScrubberError {
    /// Represents a failure to execute Kraken2
    #[error("execution error - is Kraken2 installed?")]
    KrakenError,
    /// Represents a failure in the correct length of the input file vector
    #[error("input file error - incorrect number of input files")]
    FileNumberError,
    /// Represents a failure to convert a PathBuf converted OsStr into a String
    /// This is because into_string() returns Result<String, OsString)
    #[error("input file error - incorrect format of the input file path: are there non-standard characters?")]
    InvalidFilePathConversion,
    /// Represents all other cases of `std::io::Error`.
    #[error(transparent)]
    IOError(#[from] std::io::Error),
    /// Indicates that the working directory already exists
    #[error("working directory exists ({0})")]
    WorkdirExists(String),
    /// Indicates a failure to create the working directory
    #[error("working directory path could not be created ({0})")]
    WorkdirCreateFailure(String),
    /// Indicates a failure to obtain an absolute path
    #[error("absolute path could not be obtained ({0})")]
    AbsolutePathFailure(String)
}

pub struct Scrubber {
    workdir: PathBuf
}

impl Scrubber {
    pub fn new(workdir: Option<PathBuf>) -> Result<Self, ScrubberError> {
        let _workdir = check_or_create_workdir(workdir)?;
        Ok(Self { workdir: _workdir })
    }
    pub fn run_kraken(
        &self,
        input: &Vec<PathBuf>,
        kraken_db: PathBuf,
        kraken_threads: u32,
    ) -> Result<(), ScrubberError>{
        
        // Safely build the arguments for Kraken2

        let kraken_db_path = kraken_db.into_os_string().into_string().map_err(|_| ScrubberError::InvalidFilePathConversion)?;

        let file_args = match input.len() {
            2 => {
                let file1 = input[0].clone().into_os_string().into_string().map_err(|_| ScrubberError::InvalidFilePathConversion)?;
                let file2 = input[1].clone().into_os_string().into_string().map_err(|_| ScrubberError::InvalidFilePathConversion)?;
                [Some(file1), Some(file2)]
            },
            1 => [Some(input[0].clone().into_os_string().into_string().map_err(|_| ScrubberError::InvalidFilePathConversion)?), None],
            _ => return Err(ScrubberError::FileNumberError),
        };

        let paired_arg = match input.len() {
            2 => Some("--paired"),
            1 => None,
            _ => return Err(ScrubberError::FileNumberError),
        };

        let kraken_threads_arg = kraken_threads.to_string();
        let mut kraken_args = Vec::from(["--threads".to_string(), kraken_threads_arg, "--db".to_string(), kraken_db_path]);

        match paired_arg {
            Some(value) => kraken_args.push(value.to_string()),
            None => {}
        };

        for file in file_args {
            match file {
                Some(value) => {
                    kraken_args.push(value)
                }
                None => {}
            }
        };

        // Run the Kraken command
        
        let output = Command::new("kraken2")
                            .args(kraken_args)
                            .current_dir(&self.workdir)
                            .output()
                            .map_err(|_| ScrubberError::KrakenError)?;

        println!("status: {}", output.status);
        println!("stdout: {}", String::from_utf8_lossy(&output.stdout));
        println!("stderr: {}", String::from_utf8_lossy(&output.stderr));

        Ok(())
    }

}

/// Checks if work directory exists and otherwise creates it - return the actual path for the application
/// 
/// # Errors
/// A [`ScrubberError::WorkdirExists`](#clierror) is returned if the directory already exists
/// A [`ScrubberError::WorkdirCreateFailure`](#clierror) is returned if the directory cannot be created
pub fn check_or_create_workdir(workdir: Option<PathBuf>) -> Result<PathBuf, ScrubberError> {
    let _workdir = match workdir {
        Some(path) => path,
        None => PathBuf::from(format!("scrubby_{:}", Local::now().format("%Y%m%dT%H%M%S")))
    };
    let abs_workdir = canonicalize(&_workdir).map_err(|_| ScrubberError::AbsolutePathFailure(format!("{:?}", _workdir)))?;
    let abs_workdir_msg = format!("{:?}", abs_workdir);
    if !&_workdir.exists(){
        create_dir_all(&_workdir).map_err(|_| ScrubberError::WorkdirCreateFailure(abs_workdir_msg))?;
        Ok(abs_workdir)
    } else {
        Err(ScrubberError::WorkdirExists(abs_workdir_msg))
    }
}