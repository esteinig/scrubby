use std::fs::File;
use std::io::BufWriter;
use std::process::{Command, Output, Stdio};
use needletail::{parse_fastx_file, FastxReader};
use niffler::get_writer;
use tempfile::{Builder, TempDir};
use std::collections::HashSet;
use std::path::PathBuf;
use rayon::iter::ParallelIterator;
use rayon::iter::IntoParallelRefIterator;
use std::ffi::OsStr;
use std::path::Path;

use crate::error::ScrubbyError;
use crate::scrubby::{Aligner, Classifier, Scrubby};
use crate::classifier::{get_taxid_reads_kraken, get_taxid_reads_metabuli, get_taxids_from_report};

pub struct SamtoolsConfig {
    filter: String,
    fastq: String,
}
impl SamtoolsConfig {
    pub fn from_scrubby(scrubby: &Scrubby) -> Self {

        let threads = scrubby.config.samtools_threads.unwrap_or(4);

        let filter = if scrubby.reverse { 
            "samtools view -hF 12 -".to_string() 
        } else { 
            "samtools view -f 12 -".to_string() 
        };

        let fastq = if scrubby.config.paired_end {
            format!(
                "samtools fastq --threads {threads} {} -c 6 -n -1 '{}' -2 '{}'", 
                match scrubby.config.unpaired { true => "", false => "-s /dev/null" }, 
                scrubby.output[0].display(), 
                scrubby.output[1].display()
            )
        } else {
            format!("samtools fastq --threads {threads} -c 6 -n -0 '{}'", scrubby.output[0].display())
        };

        Self { filter, fastq }
    }
    
    pub fn get_pipeline(&self) -> String {
        format!("{} | {}", self.filter, self.fastq)
    }
}


pub struct Cleaner {
    scrubby: Scrubby,
    samtools: SamtoolsConfig
}

impl Cleaner {
    pub fn from_scrubby(scrubby: &Scrubby) -> Result<Self, ScrubbyError> {

        let pipeline = Cleaner { 
            scrubby: scrubby.clone(), 
            samtools: SamtoolsConfig::from_scrubby(&scrubby)
        };

        if let Some(aligner) = &pipeline.scrubby.config.aligner {
            pipeline.check_aligner_dependency(aligner)?;
        } else if let Some(classifier) = &pipeline.scrubby.config.classifier {
            pipeline.check_classifier_dependency(classifier)?;
        }

        Ok(pipeline)

    }
    pub fn run_aligner(&self) -> Result<(), ScrubbyError> {
        match self.scrubby.config.aligner {
            Some(Aligner::Minimap2) => self.run_minimap2()?,
            Some(Aligner::Bowtie2) => self.run_bowtie2()?,
            Some(Aligner::Strobealign) => self.run_strobealign()?,
            None => return Err(ScrubbyError::MissingAligner),
        }
        Ok(())
    }
    pub fn run_classifier(&self) -> Result<(), ScrubbyError> {
        match self.scrubby.config.classifier {
            Some(Classifier::Kraken2) => self.run_kraken()?,
            Some(Classifier::Metabuli) => self.run_metabuli()?,
            None => return Err(ScrubbyError::MissingClassifier),
        }
        Ok(())
    }
    pub fn clean_reads(&self, read_ids: &HashSet<String>) -> Result<(), ScrubbyError> {
        
        if self.scrubby.config.paired_end {

            rayon::ThreadPoolBuilder::new().num_threads(if self.scrubby.config.needletail_parallel { 2 } else { 1 }).build()?.install(|| -> Result<(), ScrubbyError> {
                [0 as usize, 1 as usize].par_iter().map(|i| -> Result<(), ScrubbyError> {

                    let fastq_cleaner = FastqCleaner::from(
                        &self.scrubby.input[*i],
                        &self.scrubby.output[*i]
                    );
                    fastq_cleaner.clean_reads(&read_ids, self.scrubby.reverse)?;

                    Ok(())
                }).collect::<Result<Vec<_>, ScrubbyError>>()?;

                Ok(())
            })?;
        } else {
            let fastq_cleaner = FastqCleaner::from(
                &self.scrubby.input[0],
                &self.scrubby.output[0]
            );
            fastq_cleaner.clean_reads(&read_ids, self.scrubby.reverse)?;
        }
        Ok(())
    }
    fn check_aligner_dependency(&self, aligner: &Aligner) -> Result<(), ScrubbyError> {
        let command = match aligner {
            Aligner::Minimap2 => "minimap2 --version",
            Aligner::Bowtie2 => "bowtie2 --version",
            Aligner::Strobealign => "strobealign --version",
        };
        self.run_version_command(command).map_err(|_| ScrubbyError::AlignerDependencyMissing(aligner.clone()))?;

        Ok(())
    }
    fn check_classifier_dependency(&self, classifier: &Classifier) -> Result<(), ScrubbyError> {
        let command = match classifier {
            Classifier::Kraken2 => "kraken2 --version",
            Classifier::Metabuli => "metabuli",
        };
        self.run_version_command(command).map_err(|_| ScrubbyError::ClassifierDependencyMissing(classifier.clone()))?;

        Ok(())
    }
    fn run_version_command(&self, command: &str) -> Result<Output, ScrubbyError> {
        let output = Command::new("sh")
            .arg("-c")
            .arg(command)
            .output()
            .map_err(|e| ScrubbyError::CommandExecutionFailed(command.to_string(), e.to_string()))?;

        if !output.status.success() {
            return Err(ScrubbyError::CommandFailed(command.to_string(), output.status.code().unwrap_or(-1)));
        }

        Ok(output)
    }
    fn run_kraken(&self) -> Result<(), ScrubbyError> {

        let classifier_args = self.scrubby.config.classifier_args.as_deref().unwrap_or("");
        let classifier_index = self.scrubby.config.classifier_index.as_ref().ok_or(ScrubbyError::MissingClassifierIndex)?;
        
        let temp_dir = match &self.scrubby.workdir {
            Some(path) => Builder::new().tempdir_in(path)?,
            None => TempDir::new()?,
        };
        
        let kraken_reads = temp_dir.path().join("kraken.reads");
        let kraken_report = temp_dir.path().join("kraken.report");

        let cmd = if self.scrubby.config.paired_end {
            format!(
                "kraken2 --threads {} --db {} {} --paired {} {} --output {} --report {}",
                self.scrubby.threads,
                classifier_index.display(),
                classifier_args,
                self.scrubby.input[0].display(),
                self.scrubby.input[1].display(),
                kraken_reads.display(),
                kraken_report.display(),
            )
        } else {
            format!(
                "kraken2 --threads {} --db {} {} --single {} --output {} --report {}",
                self.scrubby.threads,
                classifier_index.display(),
                classifier_args,
                self.scrubby.input[0].display(),
                kraken_reads.display(),
                kraken_report.display(),
            )
        };

        self.run_command(&cmd)?;

        self.clean_reads(
            &self.parse_classifier_output(
                &kraken_report, 
                &kraken_reads
            )?
        )?;

        temp_dir.close()?;

        Ok(())
    }
    fn run_metabuli(&self) -> Result<(), ScrubbyError> {
        let classifier_args = self.scrubby.config.classifier_args.as_deref().unwrap_or("");
        let classifier_index = self.scrubby.config.classifier_index.as_ref().ok_or(ScrubbyError::MissingClassifierIndex)?;

        let temp_dir = match &self.scrubby.workdir {
            Some(path) => Builder::new().tempdir_in(path)?,
            None => TempDir::new()?,
        };

        let cmd = if self.scrubby.config.paired_end {
            format!(
                "metabuli classify --seq-mode 2 --threads {} {} {} {} {} {} {}",
                self.scrubby.threads,
                classifier_args,
                self.scrubby.input[0].display(),
                self.scrubby.input[1].display(),
                classifier_index.display(),
                temp_dir.path().display(),
                "metabuli".to_string()
            )
        } else {
            format!(
                "metabuli classify --seq-mode 3 --threads {} {} {} {} {} {}",
                self.scrubby.threads,
                classifier_args,
                self.scrubby.input[0].display(),
                classifier_index.display(),
                temp_dir.path().display(),
                "metabuli".to_string()
            )
        };

        self.run_command(&cmd)?;

        self.clean_reads(
            &self.parse_classifier_output(
                &temp_dir.path().join("metabuli_report.tsv"), 
                &temp_dir.path().join("metabuli_classifications.tsv")
            )?
        )?;

        temp_dir.close()?;

        Ok(())
    }
    fn parse_classifier_output(&self, report: &PathBuf, reads: &PathBuf) -> Result<HashSet<String>, ScrubbyError> {
        let taxids = get_taxids_from_report(
            report, 
            &self.scrubby.config.taxa, 
            &self.scrubby.config.taxa_direct
        )?;
        
        match &self.scrubby.config.classifier {
            Some(Classifier::Kraken2) => Ok(get_taxid_reads_kraken(taxids, reads)?),
            Some(Classifier::Metabuli) => Ok(get_taxid_reads_metabuli(taxids, reads)?),
            None => Err(ScrubbyError::MissingClassifier),
        }
    }
    fn run_minimap2(&self) -> Result<(), ScrubbyError> {
        let aligner_args = self.scrubby.config.aligner_args.as_deref().unwrap_or("");
        let alignment_index = self.scrubby.config.aligner_index.as_ref().ok_or(ScrubbyError::MissingAlignmentIndex)?;

        let cmd = if self.scrubby.config.paired_end {
            format!(
                "minimap2 -ax sr -m 40 --secondary=no -t {} {} '{}' '{}' '{}' | {}",
                self.scrubby.threads, 
                aligner_args, 
                alignment_index.display(), 
                self.scrubby.input[0].display(), 
                self.scrubby.input[1].display(), 
                self.samtools.get_pipeline()
            )
        } else {
            format!(
                "minimap2 -ax map-ont -m 40 --secondary=no -t {} {} '{}' '{}' | {}",
                self.scrubby.threads, 
                aligner_args, 
                alignment_index.display(), 
                self.scrubby.input[0].display(), 
                self.samtools.get_pipeline()
            )
        };
        self.run_command(&cmd)
    }
    fn run_bowtie2(&self) -> Result<(), ScrubbyError> {
        let aligner_args = self.scrubby.config.aligner_args.as_deref().unwrap_or("");
        let alignment_index = self.scrubby.config.aligner_index.as_ref().ok_or(ScrubbyError::MissingAlignmentIndex)?;
        
        let cmd = if self.scrubby.config.paired_end {
            format!(
                "bowtie2 -x '{}' -1 '{}' -2 '{}' -k 1 --mm -p {} {} | {}",
                alignment_index.display(), 
                self.scrubby.input[0].display(), 
                self.scrubby.input[1].display(), 
                self.scrubby.threads, 
                aligner_args, 
                self.samtools.get_pipeline()
            )
        } else {
            format!(
                "bowtie2 -x '{}' -U '{}' -k 1 --mm -p {} {} | {} ",
                alignment_index.display(), 
                self.scrubby.input[0].display(), 
                self.scrubby.threads, 
                aligner_args,
                self.samtools.get_pipeline()
            )
        };
        self.run_command(&cmd)
    }
    fn run_strobealign(&self) -> Result<(), ScrubbyError> {
        let aligner_args = self.scrubby.config.aligner_args.as_deref().unwrap_or("");
        let alignment_index = self.scrubby.config.aligner_index.as_ref().ok_or(ScrubbyError::MissingAlignmentIndex)?;
        
        let cmd = if self.scrubby.config.paired_end {
            format!(
                "strobealign -t {} {} '{}' '{}' '{}' | {}",
                self.scrubby.threads, 
                aligner_args, 
                alignment_index.display(), 
                self.scrubby.input[0].display(), 
                self.scrubby.input[1].display(), 
                self.samtools.get_pipeline()
            )
        } else {
            format!(
                "strobealign -t {} {} '{}' '{}' | {}",
                self.scrubby.threads, 
                aligner_args, 
                alignment_index.display(), 
                self.scrubby.input[0].display(), 
                self.samtools.get_pipeline(),
            )
        };
        self.run_command(&cmd)
    }
    fn run_command(&self, cmd: &str) -> Result<(), ScrubbyError> {
        
        log::debug!("Running command: {}", cmd);

        let status = Command::new("sh")
            .arg("-c")
            .arg(cmd)
            .stderr(Stdio::null())
            .status()
            .map_err(|e| ScrubbyError::CommandExecutionFailed(cmd.to_string(), e.to_string()))?;

        if !status.success() {
            return Err(ScrubbyError::CommandFailed(cmd.to_string(), status.code().unwrap_or(-1)));
        }

        Ok(())
    }
}


pub struct FastqCleaner {
    input: PathBuf,
    output: PathBuf
}
impl FastqCleaner {
    pub fn from(input: &PathBuf, output: &PathBuf) -> Self {
        Self { input: input.to_owned(), output: output.to_owned() }
    }
    pub fn clean_reads(&self, read_ids: &HashSet<String>, reverse: bool) -> Result<(), ScrubbyError> {

        let (mut reader, mut writer) = get_niffler_fastx_reader_writer(
            &self.input, 
            &self.output, 
            niffler::compression::Level::Six, 
            None
        )?;

        while let Some(rec) = reader.next() {
            let record = rec?;
            let id = get_id(record.id())?;

            // Depletion 
            if !reverse && !read_ids.contains(&id) {
                record.write(&mut writer, None)?;
            }
            // Extraction 
            if reverse && read_ids.contains(&id) {
                record.write(&mut writer, None)?;
            }
        };

        Ok(())
    }
}


// Constants
pub trait CompressionExt {
    fn from_path<S: AsRef<OsStr> + ?Sized>(p: &S) -> Self;
}

/// Attempts to infer the compression type from the file extension.
/// If the extension is not known, then Uncompressed is returned.
impl CompressionExt for niffler::compression::Format {
    fn from_path<S: AsRef<OsStr> + ?Sized>(p: &S) -> Self {
        let path = Path::new(p);
        match path.extension().map(|s| s.to_str()) {
            Some(Some("gz")) => Self::Gzip,
            Some(Some("bz") | Some("bz2")) => Self::Bzip,
            Some(Some("lzma") | Some("xz")) => Self::Lzma,
            _ => Self::No,
        }
    }
}

// Utility function to get a Needletail reader and Niffler compressed/uncompressed writer
pub fn get_niffler_fastx_reader_writer(
    input: &PathBuf,
    output: &PathBuf,
    compression_level: niffler::compression::Level,
    output_format: Option<niffler::compression::Format>,
) -> Result<(Box<dyn FastxReader>, Box<dyn std::io::Write>), ScrubbyError> {

    // Input output of read files includes compression detection
    let reader = parse_fastx_file(input)?;

    let file: File = File::create(output)?;
    let file_handle = BufWriter::new(file);
    let format = match output_format {
        None => niffler::Format::from_path(output),
        Some(format) => format,
    };
    let writer = get_writer(
        Box::new(file_handle), 
        format, 
        compression_level
    )?;

    Ok((reader, writer))
        
}


pub fn get_id(id: &[u8]) -> Result<String, ScrubbyError> {
    let header = std::str::from_utf8(id)?;
    let header_components = header
        .split_whitespace()
        .collect::<Vec<&str>>();
    
    if header_components.len() < 1 {
        return Err(ScrubbyError::NeedletailFastqHeader)
    }
    let id = header_components[0].to_string();

    Ok(id)
}