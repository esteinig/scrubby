use std::path::PathBuf;
use std::io::Write;
use chrono::{SecondsFormat, Local};
use clap::crate_version;
use serde::{Deserialize, Serialize};
use crate::{error::ScrubbyError, scrubby::{Aligner, Classifier, Scrubby}, utils::ReadDifference};


#[derive(Serialize, Deserialize)]
pub struct ScrubbyReport {
    pub version: String,
    pub date: String,
    pub command: String,
    pub input: Vec<PathBuf>,
    pub output: Vec<PathBuf>,
    pub total: u64,
    pub depleted: u64,
    pub extracted: u64,
    pub settings: ScrubbySettings
}
impl ScrubbyReport {
    pub fn create(scrubby: &Scrubby, json: &Option<PathBuf>, read_ids: &Option<PathBuf>, header: bool) -> Result<Self, ScrubbyError> {

        let diff = ReadDifference::from(&scrubby.input, &scrubby.output)?;

        let report = Self {
            version: crate_version!().to_string(),
            date: Local::now().to_rfc3339_opts(SecondsFormat::Millis, false),
            command: match scrubby.config.command { 
                Some(ref cmd) => cmd.to_string(), 
                None => String::new() 
            },
            input: scrubby.input.clone(),
            output: scrubby.output.clone(),
            total: diff.result.input,
            depleted: if scrubby.reverse { 0 } else { diff.result.difference },
            extracted: if scrubby.reverse { diff.result.difference } else { 0 },
            settings: ScrubbySettings::from_scrubby(&scrubby)
        };


        if let Some(read_ids) = read_ids {
            diff.write_reads(read_ids, header)?;
        }

        if let Some(json) = json {
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
    pub aligner_index: Option<PathBuf>,
    pub alignment: Option<PathBuf>,
    pub classifier_index: Option<PathBuf>,
    pub reads: Option<PathBuf>,
    pub report: Option<PathBuf>,
    pub taxa: Vec<String>,
    pub taxa_direct: Vec<String>,
    pub classifier_args: Option<String>,
    pub aligner_args: Option<String>,
    pub min_len: u64,
    pub min_cov: f64,
    pub min_mapq: u8,
    pub reverse: bool
}
impl ScrubbySettings {
    pub fn from_scrubby(scrubby: &Scrubby) -> Self {
        Self {
            aligner: scrubby.config.aligner.clone(),
            classifier: scrubby.config.classifier.clone(),
            aligner_index: scrubby.config.aligner_index.clone(),
            alignment: scrubby.config.alignment.clone(),
            classifier_index: scrubby.config.classifier_index.clone(),
            reads: scrubby.config.classifier_reads.clone(),
            report: scrubby.config.classifier_report.clone(),
            taxa: scrubby.config.taxa.clone(),
            taxa_direct: scrubby.config.taxa_direct.clone(),
            classifier_args: scrubby.config.classifier_args.clone(),
            aligner_args: scrubby.config.aligner_args.clone(),
            min_len: scrubby.config.min_query_length,
            min_cov: scrubby.config.min_query_coverage,
            min_mapq: scrubby.config.min_mapq,
            reverse: scrubby.reverse

        }
    }
}