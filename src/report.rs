use std::path::PathBuf;
use std::io::Write;
use chrono::{SecondsFormat, Utc};
use clap::crate_version;
use serde::{Deserialize, Serialize};
use crate::{error::ScrubbyError, scrubby::{Aligner, Classifier, Preset, Scrubby}, utils::ReadDifference};


#[derive(Serialize, Deserialize)]
pub struct ScrubbyReport {
    pub version: String,
    pub date: String,
    pub command: String,
    pub input: Vec<PathBuf>,
    pub output: Vec<PathBuf>,
    pub reads_in: u64,
    pub reads_out: u64,
    pub reads_removed: u64,
    pub reads_extracted: u64,
    pub settings: ScrubbySettings
}
impl ScrubbyReport {
    pub fn create(scrubby: &Scrubby, header: bool) -> Result<Self, ScrubbyError> {

        let diff = ReadDifference::new(
            &scrubby.input, 
            &scrubby.output, 
            None, 
            None
        ).compute()?;

        let report = Self {
            version: crate_version!().to_string(),
            date: Utc::now().to_rfc3339_opts(SecondsFormat::Secs, true),
            command: match scrubby.config.command { 
                Some(ref cmd) => cmd.to_string(), 
                None => String::new() 
            },
            input: scrubby.input.clone(),
            output: scrubby.output.clone(),
            reads_in: diff.reads_in,
            reads_out: diff.reads_out,
            reads_removed: if scrubby.extract { 0 } else { diff.difference },
            reads_extracted: if scrubby.extract { diff.difference } else { 0 },
            settings: ScrubbySettings::from_scrubby(&scrubby)
        };


        if let Some(read_ids) = &scrubby.read_ids {
            diff.write_read_ids( read_ids, header)?;
        }
        if let Some(json) = &scrubby.json {
            Self::to_json(&report, json)?;
        }

        Ok(report)
    }
    pub fn to_json(report: &Self, output: &PathBuf) -> Result<(), ScrubbyError> {
        let mut file = std::fs::File::create(output)?;
        let json_string = serde_json::to_string_pretty(report)?;
        file.write_all(json_string.as_bytes())?;
        Ok(())
    }
}

#[derive(Serialize, Deserialize)]
pub struct ScrubbySettings {
    pub aligner: Option<Aligner>,
    pub classifier: Option<Classifier>,
    pub index: Option<PathBuf>,
    pub alignment: Option<PathBuf>,
    pub reads: Option<PathBuf>,
    pub report: Option<PathBuf>,
    pub taxa: Vec<String>,
    pub taxa_direct: Vec<String>,
    pub classifier_args: Option<String>,
    pub aligner_args: Option<String>,
    pub preset: Option<Preset>,
    pub min_len: u64,
    pub min_cov: f64,
    pub min_mapq: u8,
    pub extract: bool
}
impl ScrubbySettings {
    pub fn from_scrubby(scrubby: &Scrubby) -> Self {
        Self {
            aligner: scrubby.config.aligner.clone(),
            classifier: scrubby.config.classifier.clone(),
            index: scrubby.config.index.clone(),
            alignment: scrubby.config.alignment.clone(),
            reads: scrubby.config.reads.clone(),
            report: scrubby.config.report.clone(),
            taxa: scrubby.config.taxa.clone(),
            taxa_direct: scrubby.config.taxa_direct.clone(),
            classifier_args: scrubby.config.classifier_args.clone(),
            aligner_args: scrubby.config.aligner_args.clone(),
            preset: scrubby.config.preset.clone(),
            min_len: scrubby.config.min_query_length,
            min_cov: scrubby.config.min_query_coverage,
            min_mapq: scrubby.config.min_mapq,
            extract: scrubby.extract
        }
    }
}