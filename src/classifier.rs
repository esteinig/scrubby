//! This module provides functionalities for parsing and processing taxonomic data 
//! from Kraken and Metabuli output files. It includes structures and functions 
//! for handling taxonomic levels, counting reads, and extracting taxonomic identifiers.

use anyhow::Result;
use std::collections::HashMap;
use std::collections::HashSet;
use std::fmt;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::PathBuf;

use crate::error::ScrubbyError;

/// Enumeration representing taxonomic levels.
///
/// Each variant is associated with an integer value for comparison of levels.
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
pub enum TaxonomicLevel {
    None,
    Unclassified,
    NoRank,
    Root,
    Domain,
    Kingdom,
    Phylum,
    Class,
    Order,
    Family,
    Genus,
    Species,
    Unspecified,
}

impl fmt::Display for TaxonomicLevel {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            TaxonomicLevel::None => write!(f, "None"),
            TaxonomicLevel::Unclassified => write!(f, "Unclassified"),
            TaxonomicLevel::NoRank => write!(f, "NoRank"),
            TaxonomicLevel::Root => write!(f, "Root"),
            TaxonomicLevel::Domain => write!(f, "Domain"),
            TaxonomicLevel::Kingdom => write!(f, "Kingdom"),
            TaxonomicLevel::Phylum => write!(f, "Phylum"),
            TaxonomicLevel::Class => write!(f, "Class"),
            TaxonomicLevel::Order => write!(f, "Order"),
            TaxonomicLevel::Family => write!(f, "Family"),
            TaxonomicLevel::Genus => write!(f, "Genus"),
            TaxonomicLevel::Species => write!(f, "Species"),
            TaxonomicLevel::Unspecified => write!(f, "Unspecified"),
        }
    }
}

/// Structure for counting reads associated with taxa.
///
/// Provides implementation for updating and formatting the counts.
#[derive(Debug, Clone)]
pub struct TaxonCounts {
    taxa: HashMap<String, HashMap<String, u64>>,
}

impl TaxonCounts {
    /// Creates a new `TaxonCounts` instance.
    ///
    /// # Example
    ///
    /// ```
    /// use scrubby::TaxonCounts;
    /// let taxon_counts = TaxonCounts::new();
    /// ```
    pub fn new() -> Self {
        TaxonCounts {
            taxa: HashMap::new(),
        }
    }

    /// Updates the count of reads for a given taxon.
    ///
    /// # Arguments
    ///
    /// * `tax_name` - The name of the taxon.
    /// * `tax_parent` - The parent taxon.
    /// * `tax_reads` - The number of reads to add.
    ///
    /// # Example
    ///
    /// ```
    /// taxon_counts.update("GenusA".to_string(), "FamilyA".to_string(), 10);
    /// ```
    pub fn update(&mut self, tax_name: String, tax_parent: String, tax_reads: u64) {
        self.taxa
            .entry(tax_parent)
            .and_modify(|sub_counts| {
                sub_counts
                    .entry(tax_name.clone())
                    .and_modify(|tax_count| *tax_count += tax_reads)
                    .or_insert(tax_reads);
            })
            .or_insert(HashMap::from([(tax_name.clone(), tax_reads)]));
    }
}

/// Parses the Kraken output report file to extract taxonomic identifiers.
///
/// This function implements the logic to parse the report file, extracting specified taxa
/// by name or identifier, including sub-levels where appropriate.
///
/// # Arguments
///
/// * `kraken_report` - The path to the Kraken taxonomic report file.
/// * `kraken_taxa` - A list of taxa names or identifiers to extract.
/// * `kraken_taxa_direct` - A list of taxa names or identifiers to extract directly.
///
/// # Returns
///
/// * `Result<HashSet<String>, ScrubbyError>` - A set of extracted taxonomic identifiers.
///
/// # Example
///
/// ```
/// let taxids = get_taxids_from_report(&report_path, &vec!["Eukaryota".to_string()], &vec![]).unwrap();
/// ```
pub fn get_taxids_from_report(
    kraken_report: &PathBuf,
    kraken_taxa: &[String],
    kraken_taxa_direct: &[String],
) -> Result<HashSet<String>, ScrubbyError> {

    let report = BufReader::new(File::open(kraken_report)?);

    let kraken_taxa: Vec<String> = kraken_taxa.iter().map(|x| x.trim().to_string()).collect();
    let kraken_taxa_direct: Vec<String> = kraken_taxa_direct.iter().map(|x| x.trim().to_string()).collect();

    let mut taxids: HashSet<String> = HashSet::new();
    let mut tax_counts: TaxonCounts = TaxonCounts::new();

    let mut extract_taxlevel: TaxonomicLevel = TaxonomicLevel::None;
    let mut extract_parent: String = String::from("");

    'report: for line in report.lines() {
        let record: KrakenReportRecord = KrakenReportRecord::from_str(line?)?;
        let tax_level = get_tax_level(&record);

        if kraken_taxa_direct.contains(&record.tax_name) || kraken_taxa_direct.contains(&record.tax_id) {
            log::debug!(
                "Detected direct taxon to deplete ({} : {} : {} : {})",
                &tax_level.to_string(),
                &record.tax_level,
                &record.tax_id,
                &record.tax_name
            );
            taxids.insert(record.tax_id.clone());
            tax_counts.update(record.tax_name.clone(), record.tax_name.clone(), record.reads_direct);
        }

        if tax_level < TaxonomicLevel::Domain {
            log::debug!(
                "Found taxon above `Domain` - ignoring level and sub-levels  ({} : {} : {} : {})",
                &tax_level.to_string(),
                &record.tax_level,
                &record.tax_id,
                &record.tax_name
            );
            continue 'report;
        }

        if kraken_taxa.contains(&record.tax_name) || kraken_taxa.contains(&record.tax_id) {
            log::debug!(
                "Detected taxon level ({} : {} : {} : {})",
                &tax_level.to_string(),
                &record.tax_level,
                &record.tax_id,
                &record.tax_name
            );
            extract_taxlevel = tax_level;
            extract_parent = record.tax_name.clone();

            log::debug!(
                "Setting taxon level for parsing sub-levels to {} ({})",
                extract_taxlevel.to_string(),
                &record.tax_name
            );
            if record.reads_direct > 0 {
                taxids.insert(record.tax_id);
                tax_counts.update(record.tax_name.clone(), record.tax_name.clone(), record.reads_direct);
            }
        } else {
            if extract_taxlevel == TaxonomicLevel::None {
                log::debug!(
                    "Ignoring record ({} : {} : {} : {} : {})",
                    &tax_level.to_string(),
                    &record.tax_level,
                    &record.tax_id,
                    &record.tax_name,
                    &record.reads_direct
                );
                continue 'report;
            }
            if (tax_level <= extract_taxlevel) && (record.tax_level.len() == 1) {
                log::debug!(
                    "Detected taxon level for sub-level reset ({} : {} : {} : {})",
                    &tax_level.to_string(),
                    &record.tax_level,
                    &record.tax_id,
                    &record.tax_name
                );
                extract_taxlevel = TaxonomicLevel::None;
            } else {
                if record.reads_direct > 0 {
                    log::debug!(
                        "Detected taxon sub-level with reads ({} : {} : {} : {})",
                        &tax_level.to_string(),
                        &record.tax_level,
                        &record.tax_id,
                        &record.tax_name
                    );
                    taxids.insert(record.tax_id);
                    match extract_parent.as_str() {
                        "" => return Err(ScrubbyError::KrakenReportTaxonParent),
                        _ => tax_counts.update(record.tax_name.clone(), extract_parent.clone(), record.reads_direct),
                    }
                }
            }
        }
    }

    let num_taxids = taxids.len();
    let num_taxids_chars = num_taxids.to_string().len();

    log::debug!("{}", "=".repeat(55 + num_taxids_chars));
    log::debug!(
        "{} taxonomic levels with directly assigned reads detected",
        num_taxids
    );
    log::debug!("{}", "=".repeat(55 + num_taxids_chars));

    let mut reads = 0;
    for (parent, subtaxa) in tax_counts.taxa.iter() {
        for (child, count) in subtaxa {
            log::debug!("{} :: {} ({})", parent, child, count);
            reads += count;
        }
    }

    let num_reads_chars = reads.to_string().len();
    log::debug!("{}", "=".repeat(46 + num_reads_chars));
    log::debug!("{} directly assigned reads collected from report", reads);
    log::debug!("{}", "=".repeat(46 + num_reads_chars));

    Ok(taxids)
}

/// Extracts read identifiers for given taxonomic identifiers from a Kraken reads file.
///
/// # Arguments
///
/// * `taxids` - A set of taxonomic identifiers.
/// * `kraken_reads` - The path to the Kraken reads file.
///
/// # Returns
///
/// * `Result<HashSet<String>, ScrubbyError>` - A set of read identifiers matching the given taxonomic identifiers.
///
/// # Example
///
/// ```
/// let read_ids = get_taxid_reads_kraken(taxids, &reads_path).unwrap();
/// ```
pub fn get_taxid_reads_kraken(
    taxids: HashSet<String>,
    kraken_reads: &PathBuf,
) -> Result<HashSet<String>, ScrubbyError> {
    let mut reads: HashSet<String> = HashSet::new();

    if !kraken_reads.exists() {
        return Ok(reads);
    }

    let file = BufReader::new(File::open(&kraken_reads)?);
    for line in file.lines() {
        let record: KrakenReadRecord = KrakenReadRecord::from_str(line?)?;
        if taxids.contains(&record.tax_id) {
            reads.insert(record.read_id.clone());
        }
    }

    log::debug!("{} matching classified reads were detected", reads.len());
    Ok(reads)
}

/// Extracts read identifiers for given taxonomic identifiers from a Metabuli reads file.
///
/// # Arguments
///
/// * `taxids` - A set of taxonomic identifiers.
/// * `metabuli_reads` - The path to the Metabuli reads file.
///
/// # Returns
///
/// * `Result<HashSet<String>, ScrubbyError>` - A set of read identifiers matching the given taxonomic identifiers.
///
/// # Example
///
/// ```
/// let read_ids = get_taxid_reads_metabuli(taxids, &reads_path).unwrap();
/// ```
pub fn get_taxid_reads_metabuli(
    taxids: HashSet<String>,
    metabuli_reads: &PathBuf,
) -> Result<HashSet<String>, ScrubbyError> {
    let mut reads: HashSet<String> = HashSet::new();

    if !metabuli_reads.exists() {
        return Ok(reads);
    }

    let file = BufReader::new(File::open(&metabuli_reads)?);
    for line in file.lines() {
        let record: MetabuliReadRecord = MetabuliReadRecord::from_str(line?)?;
        if taxids.contains(&record.tax_id) {
            reads.insert(record.read_id.clone());
        }
    }

    log::debug!("{} matching classified reads were detected", reads.len());
    Ok(reads)
}

/// Utility function to extract the taxonomic level from a Kraken report record.
///
/// # Arguments
///
/// * `record` - A reference to a `KrakenReportRecord`.
///
/// # Returns
///
/// * `TaxonomicLevel` - The extracted taxonomic level.
///
/// # Example
///
/// ```
/// let tax_level = get_tax_level(&record);
/// ```
pub fn get_tax_level(record: &KrakenReportRecord) -> TaxonomicLevel {
    let tax_level_str = &record.tax_level;

    if tax_level_str.starts_with('U') {
        TaxonomicLevel::Unclassified
    } else if tax_level_str.starts_with("no rank") {
        TaxonomicLevel::NoRank
    } else if tax_level_str.starts_with('R') {
        TaxonomicLevel::Root
    } else if tax_level_str.starts_with('D') || tax_level_str.starts_with("superkingdom") {
        TaxonomicLevel::Domain
    } else if tax_level_str.starts_with('K') || tax_level_str.starts_with("kingdom") {
        TaxonomicLevel::Kingdom
    } else if tax_level_str.starts_with('P') || tax_level_str.starts_with("phylum") {
        TaxonomicLevel::Phylum
    } else if tax_level_str.starts_with('C') || tax_level_str.starts_with("class") {
        TaxonomicLevel::Class
    } else if tax_level_str.starts_with('O') || tax_level_str.starts_with("order") {
        TaxonomicLevel::Order
    } else if tax_level_str.starts_with('F') || tax_level_str.starts_with("family") {
        TaxonomicLevel::Family
    } else if tax_level_str.starts_with('G') || tax_level_str.starts_with("genus") {
        TaxonomicLevel::Genus
    } else if tax_level_str.starts_with('S') || tax_level_str.starts_with("species") {
        TaxonomicLevel::Species
    } else {
        TaxonomicLevel::Unspecified
    }
}

/// Structure representing a Kraken read classification record.
#[derive(Debug, Clone)]
pub struct KrakenReadRecord {
    pub classified: bool,
    pub read_id: String,
    pub tax_id: String,
    pub read_len: String,
    pub annotation: String,
}

impl KrakenReadRecord {
    /// Creates a `KrakenReadRecord` instance from a tab-separated string.
    ///
    /// # Arguments
    ///
    /// * `kraken_line` - A string containing the tab-separated fields of a Kraken read record.
    ///
    /// # Returns
    ///
    /// * `Result<KrakenReadRecord, ScrubbyError>` - The created `KrakenReadRecord` instance.
    ///
    /// # Example
    ///
    /// ```
    /// let record = KrakenReadRecord::from_str("C\tread1\t12345\t100\tannotation".to_string()).unwrap();
    /// ```
    pub fn from_str(kraken_line: String) -> Result<Self, ScrubbyError> {
        let fields: Vec<&str> = kraken_line.split('\t').collect();

        let _classified = match fields[0] {
            "U" => false,
            "C" => true,
            _ => false,
        };

        let record = Self {
            classified: _classified,
            read_id: fields[1].trim().to_string(),
            tax_id: fields[2].trim().to_string(),
            read_len: fields[3].trim().to_string(),
            annotation: fields[4].trim().to_string(),
        };

        Ok(record)
    }
}

/// Structure representing a Kraken report record.
#[derive(Debug, Clone)]
pub struct KrakenReportRecord {
    pub fraction: String,
    pub reads: u64,
    pub reads_direct: u64,
    pub tax_level: String,
    pub tax_id: String,
    pub tax_name: String,
}

impl KrakenReportRecord {
    /// Creates a `KrakenReportRecord` instance from a tab-separated string.
    ///
    /// # Arguments
    ///
    /// * `report_line` - A string containing the tab-separated fields of a Kraken report record.
    ///
    /// # Returns
    ///
    /// * `Result<KrakenReportRecord, ScrubbyError>` - The created `KrakenReportRecord` instance.
    ///
    /// # Example
    ///
    /// ```
    /// let record = KrakenReportRecord::from_str("0.05\t100\t50\tS\t12345\ttaxon_name".to_string()).unwrap();
    /// ```
    pub fn from_str(report_line: String) -> Result<Self, ScrubbyError> {
        let fields: Vec<&str> = report_line.split('\t').collect();

        let record = Self {
            fraction: fields[0].to_string(),
            reads: fields[1]
                .parse::<u64>()
                .map_err(|_| ScrubbyError::KrakenReportReadFieldConversion)?,
            reads_direct: fields[2]
                .parse::<u64>()
                .map_err(|_| ScrubbyError::KrakenReportDirectReadFieldConversion)?,
            tax_level: fields[3].trim().to_string(),
            tax_id: fields[4].trim().to_string(),
            tax_name: fields[5].trim().to_string(),
        };

        Ok(record)
    }
}

/// Structure representing a Metabuli read classification record.
#[derive(Debug, Clone)]
pub struct MetabuliReadRecord {
    pub classified: bool,
    pub read_id: String,
    pub tax_id: String,
    pub read_len: String,
    pub dna_score: String,
    pub rank: String,
    pub annotation: String,
}

impl MetabuliReadRecord {
    /// Creates a `MetabuliReadRecord` instance from a tab-separated string.
    ///
    /// # Arguments
    ///
    /// * `metabuli_line` - A string containing the tab-separated fields of a Metabuli read record.
    ///
    /// # Returns
    ///
    /// * `Result<MetabuliReadRecord, ScrubbyError>` - The created `MetabuliReadRecord` instance.
    ///
    /// # Example
    ///
    /// ```
    /// let record = MetabuliReadRecord::from_str("1\tread1\t12345\t100\t80.5\tgenus\tannotation".to_string()).unwrap();
    /// ```
    pub fn from_str(metabuli_line: String) -> Result<Self, ScrubbyError> {
        let fields: Vec<&str> = metabuli_line.split('\t').collect();

        let _classified = match fields[0] {
            "0" => false,
            "1" => true,
            _ => false,
        };

        let record = Self {
            classified: _classified,
            read_id: fields[1].trim().to_string(),
            tax_id: fields[2].trim().to_string(),
            read_len: fields[3].trim().to_string(),
            dna_score: fields[4].trim().to_string(),
            rank: fields[5].trim().to_string(),
            annotation: fields[6].trim().to_string(),
        };

        Ok(record)
    }
}
