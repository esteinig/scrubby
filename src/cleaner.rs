use crate::scrubby::{Aligner, Classifier, Scrubby};
use crate::error::ScrubbyError;
use std::process::{Command, Output, Stdio};

pub struct Cleaner {
    scrubby: Scrubby,
    paired_end: bool,
    samtools_filter: String,
    samtools_fastq: String
}

impl Cleaner {
    pub fn from_scrubby(scrubby: Scrubby) -> Result<Self, ScrubbyError> {

        let paired_end = scrubby.input.len() == 2;

        let samtools_filter = if scrubby.reverse { 
            "samtools view -hF 12 -".to_string() 
        } else { 
            "samtools view -f 12 -".to_string() 
        };

        let samtools_fastq = if paired_end {
            format!(
                "samtools fastq --threads 4 {} -c 6 -n -1 '{}' -2 '{}'", 
                match scrubby.config.unpaired { true => "", false => "-s /dev/null" }, 
                scrubby.output[0].display(), 
                scrubby.output[1].display()
            )
        } else {
            format!("samtools fastq --threads 4 -c 6 -n -0 '{}'", scrubby.output[0].display())
        };

        let pipeline = Cleaner { scrubby, paired_end, samtools_filter, samtools_fastq };

        if let Some(aligner) = &pipeline.scrubby.config.aligner {
            pipeline.check_aligner_dependency(aligner)?;
        } else if let Some(classifier) = &pipeline.scrubby.config.classifier {
            pipeline.check_classifier_dependency(classifier)?;
        }

        Ok(pipeline)

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
    pub fn run_alignment(&self) -> Result<(), ScrubbyError> {
        match self.scrubby.config.aligner {
            Some(Aligner::Minimap2) => self.run_minimap2()?,
            Some(Aligner::Bowtie2) => self.run_bowtie2()?,
            Some(Aligner::Strobealign) => self.run_strobealign()?,
            None => return Err(ScrubbyError::MissingAligner),
        }
        Ok(())
    }
    fn run_minimap2(&self) -> Result<(), ScrubbyError> {
        let aligner_args = self.scrubby.config.aligner_args.as_deref().unwrap_or("");
        let alignment_index = self.scrubby.config.aligner_index.as_ref().ok_or(ScrubbyError::MissingAlignmentIndex)?;

        let cmd = if self.paired_end {
            format!(
                "minimap2 -ax sr -m 40 --secondary=no -t {} {} '{}' '{}' '{}' | {} | {}",
                self.scrubby.threads, aligner_args, alignment_index.display(), self.scrubby.input[0].display(), self.scrubby.input[1].display(), self.samtools_filter, self.samtools_fastq
            )
        } else {
            format!(
                "minimap2 -ax map-ont -m 40 --secondary=no -t {} {} '{}' '{}' | {} | {}",
                self.scrubby.threads, aligner_args, alignment_index.display(), self.scrubby.input[0].display(), self.samtools_filter, self.samtools_fastq
            )
        };
        self.run_command(&cmd)
    }
    fn run_bowtie2(&self) -> Result<(), ScrubbyError> {
        let aligner_args = self.scrubby.config.aligner_args.as_deref().unwrap_or("");
        let alignment_index = self.scrubby.config.aligner_index.as_ref().ok_or(ScrubbyError::MissingAlignmentIndex)?;
        
        let cmd = if self.paired_end {
            format!(
                "bowtie2 -x '{}' -1 '{}' -2 '{}' -k 1 --mm -p {} {} | {} | {}",
                alignment_index.display(), self.scrubby.input[0].display(), self.scrubby.input[1].display(), self.scrubby.threads, aligner_args, self.samtools_filter, self.samtools_fastq
            )
        } else {
            format!(
                "bowtie2 -x '{}' -U '{}' -k 1 --mm -p {} {} | {} | {}",
                alignment_index.display(), self.scrubby.input[0].display(), self.scrubby.threads, aligner_args, self.samtools_filter, self.samtools_fastq
            )
        };
        self.run_command(&cmd)
    }
    fn run_strobealign(&self) -> Result<(), ScrubbyError> {
        let aligner_args = self.scrubby.config.aligner_args.as_deref().unwrap_or("");
        let alignment_index = self.scrubby.config.aligner_index.as_ref().ok_or(ScrubbyError::MissingAlignmentIndex)?;
        
        let cmd = if self.paired_end {
            format!(
                "strobealign -t {} {} '{}' '{}' '{}' | {} | {}",
                self.scrubby.threads, aligner_args, alignment_index.display(), self.scrubby.input[0].display(), self.scrubby.input[1].display(), self.samtools_filter, self.samtools_fastq
            )
        } else {
            format!(
                "strobealign -t {} {} '{}' '{}' | {} | {}",
                self.scrubby.threads, aligner_args, alignment_index.display(), self.scrubby.input[0].display(), self.samtools_filter, self.samtools_fastq
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
