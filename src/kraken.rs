use anyhow::Result;
use std::collections::HashMap;
use std::collections::HashSet;
use std::fmt;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};

use crate::metabuli::MetabuliReadRecord;
use crate::scrub::ScrubberError;
use crate::utils::get_file_strings_from_input;

/// Builds the Kraken2 command from the input configuration
///
/// # Errors
/// A [`ScrubberError::InvalidFilePathConversion`](#scrubbererror) is returned if one of the input paths could not be converted to a string
/// A [`ScrubberError::FileNumberError`](#scrubbererror) is returned if for some arcane reason the input file vector is not the correct length
pub fn get_kraken_command(
    input: &Vec<PathBuf>,
    db_path: &Path,
    db_name: &str,
    db_idx: &usize,
    threads: &u32,
) -> Result<Vec<String>, ScrubberError> {
    let kraken_db_path = db_path
        .to_path_buf()
        .into_os_string()
        .into_string()
        .map_err(|_| ScrubberError::InvalidFilePathConversion)?;
    let kraken_threads_arg = threads.to_string();

    let file_arg = get_file_strings_from_input(input)?;

    let paired_arg = match input.len() {
        2 => Some("--paired"),
        1 => None,
        _ => return Err(ScrubberError::FileNumberError),
    };

    let mut kraken_args = Vec::from([
        "--threads".to_string(),
        kraken_threads_arg,
        "--db".to_string(),
        kraken_db_path,
        "--output".to_string(),
        format!("{}-{}.kraken", db_idx, db_name),
        "--report".to_string(),
        format!("{}-{}.report", db_idx, db_name),
    ]);

    if let Some(value) = paired_arg {
        kraken_args.push(value.to_string())
    };

    for file in file_arg.iter().flatten() {
        kraken_args.push(file.to_owned())
    }

    Ok(kraken_args)
}

/// Taxonomic level enumeration
///
/// Provides an integer value for comparison of levels
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
pub enum TaxonomicLevel {
    None,
    Unclassified,
    Root,
    Superkingdom,
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
            TaxonomicLevel::Root => write!(f, "Root"),
            TaxonomicLevel::Superkingdom => write!(f, "Superkingdom"),
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

/// A struct to count reads for taxa provided by user
///
/// Provides implementation for formatting the  HashMaps
/// associated with the counts.
///
/// Keys for the HashMaps are always the taxonomic names, as this
/// construct is meant solely for formatting informative outputs and
/// summaries of depletion or extraction.
///   
#[derive(Debug, Clone)]
pub struct TaxonCounts {
    taxa: HashMap<String, HashMap<String, u64>>,
}
impl TaxonCounts {
    pub fn new() -> Self {
        TaxonCounts {
            taxa: HashMap::new(),
        }
    }
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

/// Parse the Kraken output report file
///
/// This functions implements the logic to parse the report file, and extract any
/// directly specified taxon by name or identifier. It allows for extraction of any
/// taxon sub-leves of taxa by name or identifier, for example when specifying the
/// domain `Eukaryota` it will parse all sub-levels of `Eukaryota` until the next full
/// domain level (`D`) or above (`R` or `U`) is reached (which prevents sub-specifications
/// of domains to be excluded, such as `D1` or `D2`)
///
/// It should be ensured that the underlying database taxonomy does not incorrectly specify
/// sub-levels as the triggering level - for example, when specifying `Eukaryota` the
/// 16S rRNA SILVA database incorrectly specifies `Holozoa` at the domain level, so it
/// should be included in the taxa to deplete to ensure all sequences below `Eukaryota`'
/// and `Holozoa` are excluded.
///
/// Only taxa with directly assigned reads are included in the final set of taxids.
///
/// # Errors
/// A [`ScrubberError::KrakenReportReadFieldConversion`](#scrubbererror) is returned if the read field in the report file cannot be converted into `u64`
/// A [`ScrubberError::KrakenReportDirectReadFieldConversion`](#scrubbererror) is returned if the direct read field in the report file cannot be converted into `u64`
pub fn get_taxids_from_report(
    // Kraken taxonomic report
    kraken_report: PathBuf,
    kraken_taxa: &[String],
    kraken_taxa_direct: &[String],
) -> Result<HashSet<String>, ScrubberError> {

    let report = BufReader::new(File::open(kraken_report)?);

    // Make sure no trailign whitespaces are input by user - these are taxon names or taxids to deplete
    let kraken_taxa: Vec<String> = kraken_taxa.iter().map(|x| x.trim().to_string()).collect();
    let kraken_taxa_direct: Vec<String> = kraken_taxa_direct
        .iter()
        .map(|x| x.trim().to_string())
        .collect();

    let mut taxids: HashSet<String> = HashSet::new();
    let mut tax_counts: TaxonCounts = TaxonCounts::new();

    let mut extract_taxlevel: TaxonomicLevel = TaxonomicLevel::None; // make sure this makes sense to initialize
    let mut extract_parent: String = String::from("");

    'report: for line in report.lines() {
        // Iterate over the lines in the report file - it is not known which domain comes
        // first, so we set a flag when we want to extract identifiers, until the next taxonomic
        // level report line of the same kind or above, when the requested domains are checked again
        let record: KrakenReportRecord = KrakenReportRecord::from_str(line?)?;

        let tax_level = get_tax_level(&record);

        // Add the direct taxon identifier if the record matches by taxonomic name or identifer
        // (and has reads assigned directly) - this is always the case when the direct taxon is found
        if kraken_taxa_direct.contains(&record.tax_name)
            || kraken_taxa_direct.contains(&record.tax_id)
        {
            log::info!(
                "Detected direct taxon to deplete ({} : {} : {} : {})",
                &tax_level.to_string(),
                &record.tax_level,
                &record.tax_id,
                &record.tax_name
            );
            taxids.insert(record.tax_id.clone());
            tax_counts.update(
                record.tax_name.clone(),
                record.tax_name.clone(),
                record.reads_direct,
            )
        }

        // If taxon level is above Domain - do not allow it to be processed with the block statement below
        // this is to prevent failure of the logic implemented to parse the report sub-levels of a given
        // taxonomic name or identifier

        if tax_level < TaxonomicLevel::Domain {
            // Unclassified, Root --> all should be given directly!
            log::warn!(
                "Found taxon above `Domain` - ignoring level and sub-levels  ({} : {} : {} : {})",
                &tax_level.to_string(),
                &record.tax_level,
                &record.tax_id,
                &record.tax_name
            );
            continue 'report;
        }

        if kraken_taxa.contains(&record.tax_name) || kraken_taxa.contains(&record.tax_id) {
            log::info!(
                "Detected taxon level ({} : {} : {} : {})",
                &tax_level.to_string(),
                &record.tax_level,
                &record.tax_id,
                &record.tax_name
            );
            // If the current record is in the vector of taxa (and following sub-taxa) to deplete, switch on the extraction flag
            // and set the current tax level as a flag for stopping the extraction in subsequent records that are below or equal
            // to this tax level
            extract_taxlevel = tax_level;
            extract_parent = record.tax_name.clone();

            log::debug!(
                "Setting taxon level for parsing sub-levels to {} ({})",
                extract_taxlevel.to_string(),
                &record.tax_name
            );
            // Skip all records that do not have reads directly assigned to it!
            if record.reads_direct > 0 {
                taxids.insert(record.tax_id);
                tax_counts.update(
                    record.tax_name.clone(),
                    record.tax_name.clone(),
                    record.reads_direct,
                )
            }
        } else {
            if extract_taxlevel == TaxonomicLevel::None {
                // guard against no taxa given on initial loop
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
            // If the current record is not the depletion list, first check if the taxon level indicates we need to stop - this
            // is the case if the tax level is the same or below the extraction flagged tax level set when a record was found to
            // start depletion
            if (tax_level <= extract_taxlevel) && (record.tax_level.len() == 1) {
                //  guard against premature sub-level reset (e.g. D1, D2, D3)
                // Unset the extraction flag and reset the taxonomic level
                // to the lowest setting (None - does not ocurr in report)
                log::debug!(
                    "Detected taxon level for sub-level reset ({} : {} : {} : {})",
                    &tax_level.to_string(),
                    &record.tax_level,
                    &record.tax_id,
                    &record.tax_name
                );
                extract_taxlevel = TaxonomicLevel::None;
            } else {
                // Otherwise the taxonomic level is below the one set in the flag and the taxon should be depleted
                // Skip all records that do not have reads directly assigned to it!
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
                        "" => return Err(ScrubberError::KrakenReportTaxonParent),
                        _ => tax_counts.update(
                            record.tax_name.clone(),
                            extract_parent.clone(),
                            record.reads_direct,
                        ),
                    }
                }
            }
        }
    }

    let num_taxids = taxids.len();
    let num_taxids_chars = num_taxids.to_string().len();

    log::info!("{}", "=".repeat(55 + num_taxids_chars));
    log::info!(
        "{} taxonomic levels with directly assigned reads detected",
        num_taxids
    );
    log::info!("{}", "=".repeat(55 + num_taxids_chars));

    let mut reads = 0;
    for (parent, subtaxa) in tax_counts.taxa.iter() {
        for (child, count) in subtaxa {
            log::info!("{} :: {} ({})", parent, child, count);
            reads += count;
        }
    }

    let num_reads_chars = reads.to_string().len();
    log::info!("{}", "=".repeat(46 + num_reads_chars));
    log::info!("{} directly assigned reads collected from report", reads);
    log::info!("{}", "=".repeat(46 + num_reads_chars));

    Ok(taxids)
}

pub fn get_taxid_reads_kraken(
    taxids: HashSet<String>,
    kraken_reads: PathBuf,
) -> Result<HashSet<String>, ScrubberError> {
    
    // HashSet of read identifiers for later depletion
    let mut reads: HashSet<String> = HashSet::new();
    
    if !kraken_reads.exists() {
        // If no reads are input (empty file) the read file may not exist
        return Ok(reads)
    }

    // Extraction of read identifiers extracted from the report or added directly above
    let file = BufReader::new(File::open(&kraken_reads)?);
    for line in file.lines() {
        let record: KrakenReadRecord = KrakenReadRecord::from_str(line?)?;
        if taxids.contains(&record.tax_id) {
            reads.insert(record.read_id.clone());
        }
    }

    log::info!("{} matching classified reads were detected", reads.len());

    Ok(reads)
}


pub fn get_taxid_reads_metabuli(
    taxids: HashSet<String>,
    metabuli_reads: PathBuf,
) -> Result<HashSet<String>, ScrubberError> {
    
    // HashSet of read identifiers for later depletion
    let mut reads: HashSet<String> = HashSet::new();
    
    if !metabuli_reads.exists() {
        // If no reads are input (empty file) the read file may not exist
        return Ok(reads)
    }

    // Extraction of read identifiers extracted from the report or added directly above
    let file = BufReader::new(File::open(&metabuli_reads)?);
    for line in file.lines() {
        let record: MetabuliReadRecord = MetabuliReadRecord::from_str(line?)?;
        if taxids.contains(&record.tax_id) {
            reads.insert(record.read_id.clone());
        }
    }

    log::info!("{} matching classified reads were detected", reads.len());

    Ok(reads)
}

/// A utility function to extract a non-specific tax level into an enumeration
pub fn get_tax_level(record: &KrakenReportRecord) -> TaxonomicLevel {
    let tax_level_str = &record.tax_level;

    if tax_level_str.starts_with('U') {
        TaxonomicLevel::Unclassified
    } else if tax_level_str.starts_with('R') {
        TaxonomicLevel::Root
    } else if tax_level_str.starts_with('D') || [&String::from("superkingdom")].contains(&tax_level_str) {  // Metabuli?
        TaxonomicLevel::Domain
    } else if tax_level_str.starts_with('K') {
        TaxonomicLevel::Kingdom
    } else if tax_level_str.starts_with('P') {
        TaxonomicLevel::Phylum
    } else if tax_level_str.starts_with('C') {
        TaxonomicLevel::Class
    } else if tax_level_str.starts_with('O') {
        TaxonomicLevel::Order
    } else if tax_level_str.starts_with('F') {
        TaxonomicLevel::Family
    } else if tax_level_str.starts_with('G') {
        TaxonomicLevel::Genus
    } else if tax_level_str.starts_with('S') {
        TaxonomicLevel::Species
    } else {
        TaxonomicLevel::Unspecified
    }
}

/*
==============
Kraken Records
==============
*/

// Kraken read classification record - we are handling taxonomic identifiers
// as strings in case they are not numeric (e.g. GTDB)
#[derive(Debug, Clone)]
pub struct KrakenReadRecord {
    pub classified: bool,
    pub read_id: String,
    pub tax_id: String,
    pub read_len: String,
    pub annotation: String,
}

impl KrakenReadRecord {
    pub fn from_str(kraken_line: String) -> Result<Self, ScrubberError> {
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

// Kraken read classification record
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
    // Create a record from a parsed line
    pub fn from_str(report_line: String) -> Result<Self, ScrubberError> {
        let fields: Vec<&str> = report_line.split('\t').collect();

        let record = Self {
            fraction: fields[0].to_string(),
            reads: fields[1]
                .parse::<u64>()
                .map_err(|_| ScrubberError::KrakenReportReadFieldConversion)?,
            reads_direct: fields[2]
                .parse::<u64>()
                .map_err(|_| ScrubberError::KrakenReportDirectReadFieldConversion)?,
            tax_level: fields[3].trim().to_string(),
            tax_id: fields[4].trim().to_string(),
            tax_name: fields[5].trim().to_string(),
        };

        Ok(record)
    }
}
