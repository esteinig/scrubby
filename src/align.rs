use anyhow::Result;
use rust_htslib::{bam, bam::record::Cigar, bam::Read};
use std::collections::HashSet;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};
use std::str::from_utf8;
use thiserror::Error;
use itertools::Itertools;
use crate::scrub::ScrubberError;
use crate::utils::get_file_strings_from_input;

#[derive(Error, Debug)]
pub enum ReadAlignmentError {
    /// Indicates failure to read file
    #[error("failed to read file")]
    FileIO(#[from] std::io::Error),
    /// Indicates failure to parse a BAM file
    #[error("failed to parse records from BAM")]
    HTSLIBError(#[from] rust_htslib::errors::Error),
    /// Indicates failure to parse a record name from BAM file
    #[error("failed to parse record name from BAM")]
    UTF8Error(#[from] std::str::Utf8Error),
    /// Indicates failure to parse a target name from BAM file
    #[error("failed to parse a valid record target name from BAM")]
    TIDError(#[from] std::num::TryFromIntError),
    /// Indicates failure to parse an u64 from PAF
    #[error("failed to parse a valid integer from PAF")]
    PafRecordIntError(#[from] std::num::ParseIntError),
    /// Indicates failure to parse an u64 from PAF
    #[error("failed to parse a valid input format")]
    InputFormatError,
}

/// Builds the minimap2 command from the input configuration
///
/// # Errors
/// A [`ScrubberError::InvalidFilePathConversion`](#scrubbererror) is returned if one of the input paths could not be converted to a string
/// A [`ScrubberError::FileNumberError`](#scrubbererror) is returned if the input file vector is not the correct length
pub fn get_minimap2_command(
    input: &Vec<PathBuf>,
    index_path: &Path,
    index_name: &str,
    index_idx: &usize,
    threads: &u32,
    preset: &String,
    args: &str
) -> Result<Vec<String>, ScrubberError> {
    let minimap_index_path = index_path
        .to_path_buf()
        .into_os_string()
        .into_string()
        .map_err(|_| ScrubberError::InvalidFilePathConversion)?;
    let minimap_threads_arg = threads.to_string();

    let file_arg = get_file_strings_from_input(input)?;

    let mut minimap_args = Vec::from([
        "-t".to_string(),
        minimap_threads_arg,
        "-c".to_string(),
        "-x".to_string(),
        preset.to_string(),
    ]);
    
    let add_args = args.split_whitespace().collect_vec();
    for arg in add_args {
        minimap_args.push(arg.to_string())
    }

    for arg in [
        "-o".to_string(),
        format!("{}-{}.paf", index_idx, index_name),
        minimap_index_path
    ] {
        minimap_args.push(arg)
    }

    for file in file_arg.iter().flatten() {
        minimap_args.push(file.to_owned())
    }

    Ok(minimap_args)
}

enum StrobealignReferenceFormat {
    Fasta,
    Index,
}

/// Builds the strobealign command from the input configuration
///
/// # Errors
/// A [`ScrubberError::InvalidFilePathConversion`](#scrubbererror) is returned if one of the input paths could not be converted to a string
/// A [`ScrubberError::FileNumberError`](#scrubbererror) is returned if the input file vector is not the correct length
pub fn get_strobealign_command(
    input: &Vec<PathBuf>,
    index_path: &Path,
    index_name: &str,
    index_idx: &usize,
    threads: &u32,
    mode: &String,
    args: &str
) -> Result<Vec<String>, ScrubberError> {

    let strobealign_index_path = index_path
        .to_path_buf()
        .into_os_string()
        .into_string()
        .map_err(|_| ScrubberError::InvalidFilePathConversion)?;

    let strobealign_threads_arg = threads.to_string();
    let file_arg = get_file_strings_from_input(input)?;

    // Check index format by extension

    let mut strobealign_args = Vec::new();

    let index_format = match index_path.extension().map(|s| s.to_str()) {
        Some(Some("sti")) => StrobealignReferenceFormat::Index,
        Some(_) => StrobealignReferenceFormat::Fasta,
        _ => return Err(ScrubberError::StrobealignReferenceExtension),
    };

    match index_format {
        StrobealignReferenceFormat::Index => strobealign_args.push("--use-index".to_string()),
        StrobealignReferenceFormat::Fasta => {}
    }

    let (mode_ext, mode_flag) = match mode.as_str() {
        "map" => ("paf", Some("-x".to_string())),
        "align" => ("sam", None),
        _ => return Err(ScrubberError::StrobealignMode(mode.to_string())),
    };

    if let Some(map_mode) = mode_flag {
        strobealign_args.push(map_mode)
    };

    let add_args = args.split_whitespace().collect_vec();
    for arg in add_args {
        strobealign_args.push(arg.to_string())
    }

    log::info!("{} {:?}", strobealign_index_path, index_path);

    let strobealign_ref = match index_format { 
        StrobealignReferenceFormat::Fasta => strobealign_index_path, 
        StrobealignReferenceFormat::Index => remove_last_two_extensions(&index_path).unwrap_or(strobealign_index_path)
    };

    strobealign_args.push("-t".to_string());
    strobealign_args.push(strobealign_threads_arg);
    strobealign_args.push("-U".to_string());
    strobealign_args.push("-o".to_string());
    strobealign_args.push(format!("{}-{}.{}", index_idx, index_name, mode_ext));
    strobealign_args.push(strobealign_ref);

    for file in file_arg.iter().flatten() {
        strobealign_args.push(file.to_owned())
    }
    
    Ok(strobealign_args)
}

fn remove_last_two_extensions(file_name: &Path) -> Option<String> {
    let path = Path::new(file_name);
    
    if let Some(stem) = path.file_stem() {
        let stem_str = stem.to_string_lossy();
        let parts: Vec<&str> = stem_str.split('.').collect();
        
        if parts.len() > 2 {
            let new_stem = parts[..parts.len() - 2].join(".");
            return Some(new_stem);
        } else {
            return Some(stem_str.to_string());
        }
    }
    None
}

/*
=================
Read alignment
=================
*/

/// Struct that reads and filters the
/// input alignment, retaining the read
/// identifiers of reads to be depleted
#[derive(Debug, Clone)]
pub struct ReadAlignment {
    // Read identifiers to be depleted
    pub reads: HashSet<String>,
}

impl ReadAlignment {
    // Parse alignment by inferring format from file extension
    pub fn from(
        // Path to alignment file [PAF or "-"]
        path: &Path,
        // Minimum query alignment length
        min_qaln_len: u64,
        // Minimum query coverage
        min_qaln_cov: f64,
        // Minimum mapping quality
        min_mapq: u8,
        // Explicit alignment format
        alignment_format: Option<String>,
    ) -> Result<Self, ReadAlignmentError> {
        match alignment_format {
            Some(format) => match format.as_str() {
                "bam" => ReadAlignment::from_bam(path, min_qaln_len, min_qaln_cov, min_mapq),
                "paf" => ReadAlignment::from_paf(path, min_qaln_len, min_qaln_cov, min_mapq),
                "txt" => ReadAlignment::from_txt(path),
                _ => Err(ReadAlignmentError::InputFormatError),
            },
            None => match path.extension().map(|s| s.to_str()) {
                Some(Some("paf")) => {
                    ReadAlignment::from_paf(path, min_qaln_len, min_qaln_cov, min_mapq)
                }
                Some(Some("txt")) => ReadAlignment::from_txt(path),
                Some(Some("bam") | Some("sam") | Some("cram")) => {
                    ReadAlignment::from_bam(path, min_qaln_len, min_qaln_cov, min_mapq)
                }
                _ => Err(ReadAlignmentError::InputFormatError),
            },
        }
    }
    // Parses read identifiers from a one-column text file
    pub fn from_txt(path: &Path) -> Result<Self, ReadAlignmentError> {
        log::info!("Parsing read identifiers from text file");

        let file = BufReader::new(File::open(path)?);
        let mut target_reads: HashSet<String> = HashSet::new();
        for line in file.lines() {
            let line = line?;
            target_reads.insert(line);
        }
        log_pass_reads(&target_reads)?;

        Ok(Self {
            reads: target_reads,
        })
    }
    // Parse alignments from file
    pub fn from_paf(
        // Path to alignment file [PAF or "-"]
        path: &Path,
        // Minimum query alignment length
        min_qaln_len: u64,
        // Minimum query coverage
        min_qaln_cov: f64,
        // Minimum mapping quality
        min_mapq: u8,
    ) -> Result<Self, ReadAlignmentError> {
        log::info!("Parsing read identifiers from alignment (PAF)");

        let reader = BufReader::new(File::open(path)?);
        let mut target_reads: HashSet<String> = HashSet::new();
        for result in reader.lines() {
            let record: PafRecord = PafRecord::from_str(result?)?;
            if (record.query_aligned_length() >= min_qaln_len
                || record.query_coverage() >= min_qaln_cov)
                && record.mapq >= min_mapq
            {
                target_reads.insert(record.qname);
            }
        }
        log_pass_reads(&target_reads)?;

        Ok(Self {
            reads: target_reads,
        })
    }
    // Parse alignments from file
    pub fn from_bam(
        // Path to alignment file [SAM/BAM/CRAM or "-"]
        path: &Path,
        // Minimum query alignment length
        min_qaln_len: u64,
        // Minimum query coverage
        min_qaln_cov: f64,
        // Minimum mapping quality
        min_mapq: u8,
    ) -> Result<Self, ReadAlignmentError> {
        log::info!("Parsing read identifiers from alignment (SAM|BAM|CRAM)");

        let mut reader = bam::Reader::from_path(path)?;
        let mut target_reads: HashSet<String> = HashSet::new();
        for result in reader.records() {
            let record = result?;
            if record.is_unmapped() {
                continue;
            }
            let bam_record = BamRecord::from(&record)?;
            if (bam_record.qalen >= min_qaln_len || bam_record.query_coverage() >= min_qaln_cov)
                && bam_record.mapq >= min_mapq
            {
                target_reads.insert(bam_record.qname);
            }
        }
        log_pass_reads(&target_reads)?;

        Ok(Self {
            reads: target_reads,
        })
    }
}

fn log_pass_reads(reads: &HashSet<String>) -> Result<(), ReadAlignmentError> {
    let num_reads = reads.len();
    let num_reads_chars = num_reads.to_string().len();
    log::info!("{}", "=".repeat(44 + num_reads_chars));
    log::info!("{} reads passing filters detected in alignment", num_reads);
    log::info!("{}", "=".repeat(44 + num_reads_chars));
    Ok(())
}

/*
=================
Alignment records
=================
*/

/// Return the query alignment length from a CIGAR string
/// as the sum of all matches (M) and insertions (I).
///
/// PAF considers insertions but not deletions when computing
/// query start and end positions from which the query alignment
/// length is then calculated in PafRecord.
fn qalen_from_cigar<'a>(cigar: impl Iterator<Item = &'a Cigar>) -> u32 {
    cigar
        .map(|x| match x {
            Cigar::Match(_) => x.len(),
            Cigar::Ins(_) => x.len(),
            _ => 0,
        })
        .sum()
}

#[derive(Debug, Clone)]
pub struct BamRecord {
    /// Query sequence name.
    pub qname: String,
    /// Query sequence length.
    pub qlen: u32,
    /// Query alignment length.
    pub qalen: u64,
    /// Mapping quality (0-255; 255 for missing).
    pub mapq: u8,
}

impl BamRecord {
    /// Create a new (reduced) BamRecord from a BAM HTS LIB record
    pub fn from(record: &bam::Record) -> Result<Self, ReadAlignmentError> {
        let qname = from_utf8(record.qname())?.to_string();
        let qlen = record.seq_len() as u32;
        let mapq = record.mapq();

        let qalen = qalen_from_cigar(record.cigar().iter());

        Ok(Self {
            qname,
            qlen,
            qalen: qalen as u64,
            mapq,
        })
    }
    /// Coverage of the aligned query sequence.
    pub fn query_coverage(&self) -> f64 {
        match self.qlen == 0 {
            true => 0f64,
            false => self.qalen as f64 / self.qlen as f64,
        }
    }
}

/// PAF record without tags
#[derive(Debug, Clone)]
pub struct PafRecord {
    /// Query sequence name.
    pub qname: String,
    /// Query sequence length.
    pub qlen: u64,
    /// Query start (0-based; BED-like; closed).
    pub qstart: usize,
    /// Query end (0-based; BED-like; open).
    pub qend: usize,
    /// Relative strand: "+" or "-".
    pub strand: String,
    /// Target sequence name.
    pub tname: String,
    /// Target sequence length.
    pub tlen: u64,
    /// Target start on original strand (0-based).
    pub tstart: usize,
    /// Target end on original strand (0-based).
    pub tend: usize,
    /// Number of matching bases in the mapping.
    pub mlen: u64,
    /// Alignment block length. Number of bases, including gaps, in the mapping.
    pub blen: u64,
    /// Mapping quality (0-255; 255 for missing).
    pub mapq: u8,
}

impl PafRecord {
    // Create a record from a parsed line
    pub fn from_str(paf: String) -> Result<Self, ReadAlignmentError> {
        let fields: Vec<&str> = paf.split('\t').collect();

        let record = Self {
            qname: fields[0].to_string(),
            qlen: fields[1].parse::<u64>()?,
            qstart: fields[2].parse::<usize>()?,
            qend: fields[3].parse::<usize>()?,
            strand: fields[4].to_string(),
            tname: fields[5].to_string(),
            tlen: fields[6].parse::<u64>()?,
            tstart: fields[7].parse::<usize>()?,
            tend: fields[8].parse::<usize>()?,
            mlen: fields[9].parse::<u64>()?,
            blen: fields[10].parse::<u64>()?,
            mapq: fields[11].parse::<u8>()?,
        };

        Ok(record)
    }
    /// Length of the aligned query sequence.
    /// This is equal to the absolute value of `PafRecord.qend` - `PafRecord.qstart`.
    pub fn query_aligned_length(&self) -> u64 {
        (self.qend - self.qstart) as u64
    }
    /// Coverage of the aligned query sequence.
    /// Proportion of the query sequence involved in the alignment.
    /// This is equal to `PafRecord.query_aligned_length` - `PafRecord.qlen`
    pub fn query_coverage(&self) -> f64 {
        match self.qlen == 0 {
            true => 0f64,
            false => self.query_aligned_length() as f64 / self.qlen as f64,
        }
    }
}
