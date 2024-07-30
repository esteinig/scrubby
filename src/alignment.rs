use rust_htslib::{bam, bam::record::Cigar, bam::Read};
use serde::{Deserialize, Serialize};
use std::collections::HashSet;
use std::io::{BufRead, BufReader};
use std::path::PathBuf;
use std::str::from_utf8;

use crate::error::ScrubbyError;


#[derive(Serialize, Deserialize, Clone, Debug, clap::ValueEnum)]
pub enum AlignmentFormat {
    Sam,
    Bam,
    Cram,
    Paf,
    Txt
}

/// Struct that reads and filters the input alignment, retaining the read
/// identifiers of reads to be depleted
#[derive(Debug, Clone)]
pub struct ReadAlignment {
    pub target_reads: HashSet<String>,
}

impl ReadAlignment {
    pub fn from(
        path: &PathBuf,
        min_qaln_len: u64,
        min_qaln_cov: f64,
        min_mapq: u8,
        alignment_format: Option<AlignmentFormat>,
    ) -> Result<Self, ScrubbyError> {
        match alignment_format {
            Some(format) => match format {
                AlignmentFormat::Sam | AlignmentFormat::Bam | AlignmentFormat::Cram  => ReadAlignment::from_bam(path, min_qaln_len, min_qaln_cov, min_mapq),
                AlignmentFormat::Paf => ReadAlignment::from_paf(path, min_qaln_len, min_qaln_cov, min_mapq),
                AlignmentFormat::Txt => ReadAlignment::from_txt(path)
            },
            None => match path.extension().map(|s| s.to_str()) {
                Some(Some("paf")) | Some(Some("paf.gz")) | Some(Some("paf.xz")) | Some(Some("paf.bz")) | Some(Some("paf.bz2")) => ReadAlignment::from_paf(path, min_qaln_len, min_qaln_cov, min_mapq),
                Some(Some("txt")) |  Some(Some("txt.gz")) | Some(Some("txt.xz")) | Some(Some("txt.bz")) | Some(Some("txt.bz2")) => ReadAlignment::from_txt(path),
                Some(Some("bam") | Some("sam") | Some("cram")) => ReadAlignment::from_bam(path, min_qaln_len, min_qaln_cov, min_mapq),
                _ => Err(ScrubbyError::AlignmentInputFormatNotRecognized),
            },
        }
    }
    // Parses read identifiers from a one-column text file
    pub fn from_txt(path: &PathBuf) -> Result<Self, ScrubbyError> {
        
        let (reader, _) = niffler::from_path(path)?;
        let reader = BufReader::new(reader);

        let mut target_reads: HashSet<String> = HashSet::new();
        for line in reader.lines() {
            let line = line?;
            target_reads.insert(line);
        }

        Ok(Self {
            target_reads,
        })
    }
    // Parse alignments from file
    pub fn from_paf(
        path: &PathBuf,
        min_qaln_len: u64,
        min_qaln_cov: f64,
        min_mapq: u8,
    ) -> Result<Self, ScrubbyError> {
        let (reader, _) = niffler::from_path(path)?;
        let reader = BufReader::new(reader);

        let mut target_reads: HashSet<String> = HashSet::new();
        for result in reader.lines() {
            let record: PafRecord = PafRecord::from_str(result?)?;
            if (record.query_aligned_length() < min_qaln_len
                || record.query_coverage() < min_qaln_cov)
                && record.mapq < min_mapq
            {
                target_reads.insert(record.qname);
            }
        }

        Ok(Self {
            target_reads,
        })
    }
    // Parse alignments from file
    pub fn from_bam(
        path: &PathBuf,
        min_qaln_len: u64,
        min_qaln_cov: f64,
        min_mapq: u8,
    ) -> Result<Self, ScrubbyError> {

        let mut reader = bam::Reader::from_path(path)?;
        let mut target_reads: HashSet<String> = HashSet::new();

        for result in reader.records() {
            let record = result?;
            if record.is_unmapped() {
                continue;
            }
            let bam_record = BamRecord::from(&record)?;
            if (bam_record.qalen < min_qaln_len || bam_record.query_coverage() < min_qaln_cov)
                && bam_record.mapq < min_mapq
            {
                target_reads.insert(bam_record.qname);
            }
        }

        Ok(Self {
            target_reads,
        })
    }
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
    pub fn from(record: &bam::Record) -> Result<Self, ScrubbyError> {
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
    pub fn from_str(paf: String) -> Result<Self, ScrubbyError> {
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
    pub fn query_aligned_length(&self) -> u64 {
        (self.qend - self.qstart) as u64
    }
    /// Coverage of the aligned query sequence.
    /// Proportion of the query sequence involved in the alignment.
    pub fn query_coverage(&self) -> f64 {
        match self.qlen == 0 {
            true => 0f64,
            false => self.query_aligned_length() as f64 / self.qlen as f64,
        }
    }
}
