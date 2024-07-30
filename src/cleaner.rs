//! This module provides functionalities for cleaning and processing FASTQ files
//! using various aligners and classifiers. It includes the core structures and 
//! implementations for executing the cleaning pipeline with the Scrubby tool.

use std::process::{Command, Output, Stdio};
use tempfile::{Builder, TempDir};
use std::collections::HashSet;
use std::path::PathBuf;
use rayon::iter::ParallelIterator;
use rayon::iter::IntoParallelRefIterator;

#[cfg(mm2)] 
use crossbeam::channel;
#[cfg(mm2)]
use rayon::iter::IntoParallelIterator;
#[cfg(mm2)]
use needletail::parse_fastx_file;
#[cfg(mm2)]
use std::sync::{Arc, Mutex};

use crate::alignment::ReadAlignment;
use crate::error::ScrubbyError;
use crate::scrubby::{Aligner, Classifier, Scrubby};
use crate::classifier::{get_taxid_reads_kraken, get_taxid_reads_metabuli, get_taxids_from_report};
use crate::utils::{get_id, get_niffler_fastx_reader_writer};

/// Configuration for Samtools commands used in the cleaning process.
pub struct SamtoolsConfig {
    filter: String,
    fastq: String,
}

impl SamtoolsConfig {
    /// Constructs a new `SamtoolsConfig` from the provided `Scrubby` instance.
    ///
    /// # Arguments
    ///
    /// * `scrubby` - A reference to the `Scrubby` instance.
    ///
    /// # Example
    ///
    /// ```
    /// use scrubby::SamtoolsConfig;
    /// let samtools_config = SamtoolsConfig::from_scrubby(&scrubby_instance);
    /// ```
    pub fn from_scrubby(scrubby: &Scrubby) -> Self {
        let threads = scrubby.config.samtools_threads.unwrap_or(4);

        let filter = if scrubby.extract { 
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

    /// Constructs the complete pipeline command string combining filter and fastq commands.
    ///
    /// # Example
    ///
    /// ```
    /// let pipeline = samtools_config.get_pipeline();
    /// ```
    pub fn get_pipeline(&self) -> String {
        format!(
            "{} | {}", 
            self.filter, 
            self.fastq
        )
    }
}

/// Core structure for cleaning and processing FASTQ files.
pub struct Cleaner {
    scrubby: Scrubby,
    samtools: SamtoolsConfig,
}

impl Cleaner {
    /// Constructs a new `Cleaner` from the provided `Scrubby` instance.
    ///
    /// # Arguments
    ///
    /// * `scrubby` - A reference to the `Scrubby` instance.
    ///
    /// # Example
    ///
    /// ```
    /// use scrubby::Cleaner;
    /// let cleaner = Cleaner::from_scrubby(&scrubby_instance).unwrap();
    /// ```
    pub fn from_scrubby(scrubby: &Scrubby) -> Result<Self, ScrubbyError> {
        let pipeline = Cleaner { 
            scrubby: scrubby.clone(), 
            samtools: SamtoolsConfig::from_scrubby(&scrubby),
        };

        if let Some(aligner) = &pipeline.scrubby.config.aligner {
            pipeline.check_aligner_dependency(aligner)?;
        } else if let Some(classifier) = &pipeline.scrubby.config.classifier {
            pipeline.check_classifier_dependency(classifier)?;
        }

        Ok(pipeline)
    }

    /// Executes the aligner process.
    ///
    /// # Returns
    ///
    /// * `Result<(), ScrubbyError>` - Ok if the aligner process completes successfully, otherwise an error.
    ///
    /// # Example
    ///
    /// ```
    /// cleaner.run_aligner().unwrap();
    /// ```
    pub fn run_aligner(&self) -> Result<(), ScrubbyError> {
        match self.scrubby.config.aligner {
            Some(Aligner::Minimap2) => self.run_minimap2()?,
            Some(Aligner::Bowtie2) => self.run_bowtie2()?,
            Some(Aligner::Strobealign) => self.run_strobealign()?,
            #[cfg(mm2)]
            Some(Aligner::Minimap2Rs) => self.run_minimap2_rs()?,
            None => return Err(ScrubbyError::MissingAligner),
        }
        Ok(())
    }

    /// Executes the classifier process.
    ///
    /// # Returns
    ///
    /// * `Result<(), ScrubbyError>` - Ok if the classifier process completes successfully, otherwise an error.
    ///
    /// # Example
    ///
    /// ```
    /// cleaner.run_classifier().unwrap();
    /// ```
    pub fn run_classifier(&self) -> Result<(), ScrubbyError> {
        match self.scrubby.config.classifier {
            Some(Classifier::Kraken2) => self.run_kraken()?,
            Some(Classifier::Metabuli) => self.run_metabuli()?,
            None => return Err(ScrubbyError::MissingClassifier),
        }
        Ok(())
    }
    /// Executes the classifier output cleaning process.
    ///
    /// # Returns
    ///
    /// * `Result<(), ScrubbyError>` - Ok if the classifier output cleaning process completes successfully, otherwise an error.
    ///
    /// # Example
    ///
    /// ```
    /// cleaner.run_classifier_output().unwrap();
    /// ```
    pub fn run_classifier_output(&self) -> Result<(), ScrubbyError> {
        match self.scrubby.config.classifier {
            Some(Classifier::Kraken2) | Some(Classifier::Metabuli) => {
                self.clean_reads(
                    &self.parse_classifier_output(
                        &self.scrubby.config.classifier_report
                            .clone()
                            .ok_or(ScrubbyError::MissingClassifierClassificationReport)?, 
                        &self.scrubby.config.classifier_reads
                            .clone()
                            .ok_or(ScrubbyError::MissingClassifierReadClassfications)?
                    )?
                )?
            },
            None => return Err(ScrubbyError::MissingClassifier),
        }
        Ok(())
    }
    /// Executes the alignment output cleaning process.
    ///
    /// # Returns
    ///
    /// * `Result<(), ScrubbyError>` - Ok if the alignment output cleaning process completes successfully, otherwise an error.
    ///
    /// # Example
    ///
    /// ```
    /// cleaner.run_aligner_output().unwrap();
    /// ```
    pub fn run_aligner_output(&self) -> Result<(), ScrubbyError> {
        
        let alignment = ReadAlignment::from(
            &self.scrubby.config.alignment.clone().ok_or(ScrubbyError::MissingAlignment)?,
            self.scrubby.config.min_query_length,
            self.scrubby.config.min_query_coverage,
            self.scrubby.config.min_mapq,
            self.scrubby.config.alignment_format.clone()
        )?;

        self.clean_reads(&alignment.target_reads)?;

        Ok(())
    }
    /// Cleans reads based on the provided read IDs.
    ///
    /// # Arguments
    ///
    /// * `read_ids` - A reference to a set of read IDs to be cleaned.
    ///
    /// # Returns
    ///
    /// * `Result<(), ScrubbyError>` - Ok if the cleaning process completes successfully, otherwise an error.
    ///
    /// # Example
    ///
    /// ```
    /// let read_ids = HashSet::new();
    /// cleaner.clean_reads(&read_ids).unwrap();
    /// ```
    pub fn clean_reads(&self, read_ids: &HashSet<String>) -> Result<(), ScrubbyError> {
        if self.scrubby.config.paired_end {
            rayon::ThreadPoolBuilder::new()
                .num_threads(if self.scrubby.config.needletail_parallel { 2 } else { 1 })
                .build()?
                .install(|| -> Result<(), ScrubbyError> {
                    [0, 1].par_iter().map(|&i| {
                        let fastq_cleaner = FastqCleaner::from(&self.scrubby.input[i], &self.scrubby.output[i]);
                        fastq_cleaner.clean_reads(&read_ids, self.scrubby.extract)?;
                        Ok(())
                    }).collect::<Result<Vec<_>, ScrubbyError>>()?;
                    Ok(())
                })?;
        } else {
            let fastq_cleaner = FastqCleaner::from(&self.scrubby.input[0], &self.scrubby.output[0]);
            fastq_cleaner.clean_reads(&read_ids, self.scrubby.extract)?;
        }
        Ok(())
    }
    fn check_aligner_dependency(&self, aligner: &Aligner) -> Result<(), ScrubbyError> {
        let command = match aligner {
            Aligner::Minimap2 => "minimap2 --version",
            Aligner::Bowtie2 => "bowtie2 --version",
            Aligner::Strobealign => "strobealign --version",
            #[cfg(mm2)]
            Aligner::Minimap2Rs => return Ok(())
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
            &self.parse_classifier_output(&kraken_report, &kraken_reads)?
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
        let taxids = get_taxids_from_report(report, &self.scrubby.config.taxa, &self.scrubby.config.taxa_direct)?;
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
                "minimap2 -ax sr --secondary=no -t {} {} '{}' '{}' '{}' | {}",
                self.scrubby.threads,
                aligner_args,
                alignment_index.display(),
                self.scrubby.input[0].display(),
                self.scrubby.input[1].display(),
                self.samtools.get_pipeline()
            )
        } else {
            format!(
                "minimap2 -ax map-ont --secondary=no -t {} {} '{}' '{}' | {}",
                self.scrubby.threads,
                aligner_args,
                alignment_index.display(),
                self.scrubby.input[0].display(),
                self.samtools.get_pipeline()
            )
        };
        self.run_command(&cmd)?;

        Ok(())
    }
    #[cfg(mm2)]
    fn run_minimap2_rs(&self) -> Result<(), ScrubbyError> {

        let (sequence_sender, sequence_receiver) = channel::unbounded();
        // let (result_sender, result_receiver) = channel::unbounded();

        let aligner = minimap2::Aligner::builder();

        let aligner = if self.scrubby.config.paired_end {
            aligner.sr()
        } else {
            aligner.map_ont()
        };
        
        let aligner = aligner
            .with_cigar()
            .with_index_threads(self.scrubby.threads)
            .with_index(
                self.scrubby.config.aligner_index.clone().ok_or(
                    ScrubbyError::MissingAlignmentIndex
                )?, 
            None
        ).map_err(|err| {
            ScrubbyError::Minimap2RustAlignerBuilderFailed(err.to_string())
        })?;

        let sequence_sender = Arc::new(Mutex::new(sequence_sender));

        rayon::scope(|s| {
            let reads_1 = self.scrubby.input[0].clone();
            let sequence_sender_clone = Arc::clone(&sequence_sender);
            s.spawn(move |_| {
                let mut reader = parse_fastx_file(reads_1).expect("Failed to open read file 1");

                while let Some(rec) = reader.next() {
                    let record = rec.expect("Failed to read record");
                    sequence_sender_clone.lock().unwrap().send((get_id(record.id()).expect("Failed to get ID"), record.seq().to_vec())).expect("Failed to send sequence 1");
                }
            });

            if self.scrubby.config.paired_end {
                let reads_2 = self.scrubby.input[1].clone();
                let sequence_sender_clone = Arc::clone(&sequence_sender);
                s.spawn(move |_| {
                    let mut reader = parse_fastx_file(reads_2).expect("Failed to open read file 2");

                    while let Some(rec) = reader.next() {
                        let record = rec.expect("Failed to read record");
                        sequence_sender_clone.lock().unwrap().send((get_id(record.id()).expect("Failed to get ID"), record.seq().to_vec())).expect("Failed to send sequence 2");
                    }
                });
            }
        });

        // Drop the sequence sender to close the channel and allow the receiver to finish
        drop(sequence_sender);

        // We might want to think about a separate thread for receiving results, but we 
        // have for now included the alignment/no alignment check in the alignment
        // threads so that minimal data is returned (sequence identifiers)

        // let writer_handle = std::thread::spawn(move || -> Result<(), ScrubbyError> {
        //     for result in result_receiver.iter() {
        //         let (read_id, mappings) = result;
        //         println!("{} {:?}", read_id, mappings);
        //     }
        //     Ok(())
        // });

        // let result_sender = Arc::new(Mutex::new(result_sender));
       
        let results = rayon::ThreadPoolBuilder::new().num_threads(self.scrubby.threads).build()?.scope(|_| -> Result<_, ScrubbyError> {
            let results = sequence_receiver
                .iter()
                .collect::<Vec<_>>()
                .into_par_iter()
                .map(|(id, sequence)| -> Result<_, ScrubbyError> {
                    let mappings = aligner.map(&sequence, false, false, None, None).map_err(|err| ScrubbyError::Minimap2RustAlignmentFailed(err.to_string()))?;
                    if mappings.len() > 0 {
                        Ok(Some(id))
                    } else {
                        Ok(None)
                    }
                })
                .collect::<Vec<_>>();
                
                // .for_each_with(result_sender.clone(), |result_sender, (read_id, sequence) | {
                //     let mappings = ;
                //     result_sender.lock().unwrap().send((read_id, mappings)).expect("Failed to send mapping results")
                // });

            Ok(results)
        })?;

        // drop(result_sender);
        // writer_handle.join().expect("Writer thread panicked")?;

        let mut read_ids = HashSet::new();
        for result in results {
            let result = result?;
            if let Some(id) = result {
                read_ids.insert(id);
            }
        }

        self.clean_reads(&read_ids)?;

        Ok(())
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
        self.run_command(&cmd)?;

        Ok(())

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
        self.run_command(&cmd)?;

        Ok(())
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

/// Structure for cleaning FASTQ files based on read IDs.
pub struct FastqCleaner {
    input: PathBuf,
    output: PathBuf,
}

impl FastqCleaner {
    /// Constructs a new `FastqCleaner` from the provided input and output paths.
    ///
    /// # Arguments
    ///
    /// * `input` - A reference to the input file path.
    /// * `output` - A reference to the output file path.
    ///
    /// # Example
    ///
    /// ```
    /// use scrubby::FastqCleaner;
    /// let cleaner = FastqCleaner::from(&input_path, &output_path);
    /// ```
    pub fn from(input: &PathBuf, output: &PathBuf) -> Self {
        Self { input: input.to_owned(), output: output.to_owned() }
    }

    /// Cleans reads from the input file and writes to the output file based on the provided read IDs.
    ///
    /// # Arguments
    ///
    /// * `read_ids` - A reference to a set of read IDs to be cleaned.
    /// * `reverse` - A boolean indicating whether to reverse the cleaning process.
    ///
    /// # Returns
    ///
    /// * `Result<(), ScrubbyError>` - Ok if the cleaning process completes successfully, otherwise an error.
    ///
    /// # Example
    ///
    /// ```
    /// let read_ids = HashSet::new();
    /// cleaner.clean_reads(&read_ids, false).unwrap();
    /// ```
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
