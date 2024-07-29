use core::fmt;
use std::fs::{create_dir_all, remove_file, File};
use std::io::{BufReader, BufWriter};
use std::path::PathBuf;
use std::time::Duration;
use reqwest::blocking::Client;
use serde::{Deserialize, Serialize};
use tar::Archive;

use crate::error::ScrubbyError;
use crate::scrubby::{Aligner, Classifier};


#[derive(Serialize, Deserialize, Clone, Debug, clap::ValueEnum)]
pub enum ScrubbyIndex {
    Chm13v2
}

impl ScrubbyIndex {
    pub fn aligner_name(&self, aligner: &Aligner) -> String {
        format!("{}.{}.tar.xz", self, aligner.short_name())
    }

    pub fn classifier_name(&self, classifier: &Classifier) -> String {
        format!("{}.{}.tar.xz", self, classifier.short_name())
    }
}

impl fmt::Display for ScrubbyIndex {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            ScrubbyIndex::Chm13v2 => write!(f, "chm13v2"),
        }
    }
}

pub struct ScrubbyDownloader {
    pub outdir: PathBuf,
    pub base_url: String,
    pub username: String,
    pub password: String,
    pub client: Client,
    pub timeout: u64,
    pub indices: Vec<ScrubbyIndex>,
    pub aligners: Vec<Aligner>,
    pub classifiers: Vec<Classifier>
}

impl ScrubbyDownloader {
    pub fn list() {
        log::info!("Available index names for download: ");
        log::info!("Human T2T Reference:  {}", ScrubbyIndex::Chm13v2);
    }
    pub fn new(outdir: PathBuf, indices: Vec<ScrubbyIndex>) -> Result<Self, ScrubbyError> {
        
        ScrubbyDownloaderBuilder::new(outdir, indices).build()
    }

    pub fn download_index(&self) -> Result<(), ScrubbyError> {

        for index in &self.indices {
            for aligner in &self.aligners {
                let file_path = self.outdir.join(index.aligner_name(&aligner));
                log::info!("Downloading alignment index to file: {}", file_path.display());
                self.download(&index.aligner_name(aligner), &file_path)?;
                log::info!("Unpacking alignment index to directory: {}", self.outdir.display());
                self.unpack(&file_path, &self.outdir)?;
                log::info!("Removing alignment index: {}", file_path.display());
                remove_file(&file_path)?;
            }
            for classifier in &self.classifiers {
                let file_path = self.outdir.join(index.classifier_name(&classifier));
                log::info!("Downloading classifier index to file: {}", file_path.display());
                self.download(&index.classifier_name(classifier), &file_path)?;
                log::info!("Unpacking classifier index to directory: {}", self.outdir.display());
                self.unpack(&file_path, &self.outdir)?;
                log::info!("Removing classifier index: {}", file_path.display());
                remove_file(&file_path)?;
            }
        }

        Ok(())
    }

    fn unpack(&self, download: &PathBuf, outdir: &PathBuf) -> Result<(), ScrubbyError> {
        let file = File::open(download)?;
        let buf_reader = BufReader::new(file);
        let (reader, _compression) = niffler::get_reader(Box::new(buf_reader))?;

        let mut archive = Archive::new(reader);

        for entry in archive.entries()? {
            let mut entry = entry?;
            entry.unpack_in(outdir)?;
        }

        Ok(())
    }

    fn download(&self, file_name: &str, path: &PathBuf) -> Result<(), ScrubbyError> {
        let url = format!("{}/{}", self.base_url, file_name);
        log::info!("URL: {url} User: {} Password: {}", self.username, self.password);

        let mut response = self.client.get(&url)
            .basic_auth(&self.username, Some(&self.password))
            .timeout(Duration::from_secs(self.timeout*60))
            .send()?;

        if !response.status().is_success() {
            return Err(ScrubbyError::DownloadFailedRequest(response.status()));
        }

        let mut writer = BufWriter::new(File::create(path)?);
        response.copy_to(&mut writer)?;

        Ok(())
    }
}

pub struct ScrubbyDownloaderBuilder {
    indices: Vec<ScrubbyIndex>,
    outdir: PathBuf,
    base_url: Option<String>,
    timeout: Option<u64>,
    username: Option<String>,
    password: Option<String>,
    aligners: Option<Vec<Aligner>>,
    classifiers: Option<Vec<Classifier>>
}

impl ScrubbyDownloaderBuilder {
    pub fn new(outdir: PathBuf, indices: Vec<ScrubbyIndex>) -> Self {
        
        Self {
            outdir,
            indices,
            base_url: None,
            username: None,
            password: None,
            aligners: None,
            classifiers: None,
            timeout: None
        }
    }
    pub fn aligner<T: Into<Option<Vec<Aligner>>>>(mut self, aligner: T) -> Self {
        self.aligners = aligner.into();
        self
    }
    pub fn classifier<T: Into<Option<Vec<Classifier>>>>(mut self, classifier: T) -> Self {
        self.classifiers = classifier.into();
        self
    }
    pub fn timeout<T: Into<Option<u64>>>(mut self, timeout: T) -> Self {
        self.timeout = timeout.into();
        self
    }
    pub fn base_url<T: Into<Option<String>>>(mut self, base_url: T) -> Self {
        self.base_url = base_url.into();
        self
    }
    pub fn username<T: Into<Option<String>>>(mut self, username: T) -> Self {
        self.username = username.into();
        self
    }
    pub fn password<T: Into<Option<String>>>(mut self, password: T) -> Self {
        self.password = password.into();
        self
    }
    pub fn build(self) -> Result<ScrubbyDownloader, ScrubbyError> {
        
        if !self.outdir.exists() || !self.outdir.is_dir() {
            create_dir_all(&self.outdir)?;
        }

        let username = self.username.unwrap_or("u416706-sub1".to_string());
        let password = self.password.unwrap_or("G8tGWjHvdhUg4NGN".to_string());
        let base_url = self.base_url.unwrap_or(format!("https://{username}.your-storagebox.de/databases"));
        let aligners = self.aligners.unwrap_or(Vec::from([Aligner::Bowtie2]));
        let classifiers = self.classifiers.unwrap_or(Vec::new());
        let timeout = self.timeout.unwrap_or(30);

        Ok(ScrubbyDownloader {
            outdir: self.outdir.to_owned(),
            base_url,
            username,
            password,
            client: Client::new(),
            timeout,
            indices: self.indices.clone(),
            aligners,
            classifiers
        })
    }
}
