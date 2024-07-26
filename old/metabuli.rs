use anyhow::Result;
use itertools::Itertools;
use std::path::{Path, PathBuf};

use crate::scrub::ScrubberError;
use crate::utils::get_file_strings_from_input;

pub enum MetabuliSeqMode {
    One = 1,
    Two = 2,
    Three = 3
}
impl MetabuliSeqMode {
    pub fn from_arg(arg: &str) -> Self {
        match arg {
            "1" => Self::One,
            "2" => Self::Two,
            "3" => Self::Three,
            _ => unimplemented!("Argument {arg} is not valid for MetabuliSeqMode")
        }
    }
}

/// Builds the Metabuli command from the input configuration
///
pub fn get_metabuli_command(
    input: &Vec<PathBuf>,
    db_path: &Path,
    db_name: &str,
    db_idx: &usize,
    threads: &u32,
    seq_mode: Option<MetabuliSeqMode>,
    args: &str
) -> Result<Vec<String>, ScrubberError> {
    let metabuli_db_path = db_path
        .to_path_buf()
        .into_os_string()
        .into_string()
        .map_err(|_| ScrubberError::InvalidFilePathConversion)?;

    let metabuli_threads_arg = threads.to_string();

    let file_arg = get_file_strings_from_input(input)?;

    let mode = match seq_mode {
        Some(mode) => mode,
        None => {
            match input.len() {
                2 => MetabuliSeqMode::Two,    // paired end 
                1 => MetabuliSeqMode::Three,  // long read
                _ => return Err(ScrubberError::FileNumberError),
            }
        }
    };

    let mode_arg = (mode as i32).to_string();
    
    let mut metabuli_args = Vec::from([
        "classify".to_string(),
        "--seq-mode".to_string(),
        mode_arg,
        "--threads".to_string(),
        metabuli_threads_arg,        
    ]);

    let add_args = args.split_whitespace().collect_vec();
    for arg in add_args {
        metabuli_args.push(arg.to_string())
    }
    
    for file in file_arg.iter().flatten() {
        metabuli_args.push(file.to_owned())
    }

    metabuli_args.push(metabuli_db_path);
    
    metabuli_args.push(
        format!("{}-{}", db_idx, db_name)
    );
    metabuli_args.push(
        format!("{}-{}", db_idx, db_name)
    );

    Ok(metabuli_args)
}

