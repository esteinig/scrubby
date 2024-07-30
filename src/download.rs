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

/// Represents different indices available for Scrubby.
#[derive(Serialize, Deserialize, Clone, Debug, clap::ValueEnum)]
pub enum ScrubbyIndex {
    Chm13v2
}

impl ScrubbyIndex {
    /// Returns the aligner name formatted for the specified index.
    ///
    /// # Arguments
    ///
    /// * `aligner` - A reference to an `Aligner`.
    ///
    /// # Example
    ///
    /// ```
    /// let index = ScrubbyIndex::Chm13v2;
    /// let aligner = Aligner::new();
    /// let name = index.aligner_name(&aligner);
    /// ```
    pub fn aligner_name(&self, aligner: &Aligner) -> String {
        format!("{}.{}.tar.xz", self, aligner.short_name())
    }
    /// Returns the classifier name formatted for the specified index.
    ///
    /// # Arguments
    ///
    /// * `classifier` - A reference to a `Classifier`.
    ///
    /// # Example
    ///
    /// ```
    /// let index = ScrubbyIndex::Chm13v2;
    /// let classifier = Classifier::new();
    /// let name = index.classifier_name(&classifier);
    /// ```
    pub fn classifier_name(&self, classifier: &Classifier) -> String {
        format!("{}.{}.tar.xz", self, classifier.short_name())
    }
}

impl fmt::Display for ScrubbyIndex {
    /// Formats the ScrubbyIndex for display.
    ///
    /// # Arguments
    ///
    /// * `f` - A mutable reference to a `fmt::Formatter`.
    ///
    /// # Example
    ///
    /// ```
    /// let index = ScrubbyIndex::Chm13v2;
    /// println!("{}", index);
    /// ```
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            ScrubbyIndex::Chm13v2 => write!(f, "chm13v2"),
        }
    }
}

/// Manages the download process for Scrubby indices.
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
    /// Creates a new instance of ScrubbyDownloader.
    ///
    /// # Arguments
    ///
    /// * `outdir` - Output directory for downloaded files.
    /// * `indices` - A list of `ScrubbyIndex` to download.
    ///
    /// # Errors
    ///
    /// Returns a `ScrubbyError` if the downloader cannot be created.
    ///
    /// # Example
    ///
    /// ```
    /// let outdir = PathBuf::from("/path/to/output");
    /// let indices = vec![ScrubbyIndex::Chm13v2];
    /// let downloader = ScrubbyDownloader::new(outdir, indices);
    /// ```
    pub fn new(outdir: PathBuf, indices: Vec<ScrubbyIndex>) -> Result<Self, ScrubbyError> {
        ScrubbyDownloaderBuilder::new(outdir, indices).build()
    }
    /// Creates a new instance of ScrubbyDownloaderBuilder.
    ///
    /// # Arguments
    ///
    /// * `outdir` - Output directory for downloaded files.
    /// * `indices` - A list of `ScrubbyIndex` to download.
    ///
    /// # Example
    ///
    /// ```
    /// let outdir = PathBuf::from("/path/to/output");
    /// let indices = vec![ScrubbyIndex::Chm13v2];
    /// let builder = ScrubbyDownloader::builder(outdir, indices);
    /// ```
    pub fn builder(outdir: PathBuf, indices: Vec<ScrubbyIndex>) -> ScrubbyDownloaderBuilder {
        ScrubbyDownloaderBuilder::new(outdir, indices)
    }
    /// Lists the available index names for download.
    ///
    /// # Example
    ///
    /// ```
    /// let downloader = ScrubbyDownloader::new(outdir, indices);
    /// downloader.list();
    /// ```
    pub fn list(&self) {

        log::info!("=====================================");
        log::info!("Reference indices for anonymous users");
        log::info!("=====================================");
        log::info!("                                     ");
        log::info!("Indices are available at: {}", self.base_url);
        log::info!("Uername '{}' and password '{}'", self.username, self.password);
        log::info!("                                     ");
        log::info!("=====================================");
        log::info!("Available index names for download   ");
        log::info!("=====================================");
        log::info!("                                     ");
        log::info!("{:<16} Human T2T Reference (CHM13v2)", ScrubbyIndex::Chm13v2);
    }
    /// Downloads the specified indices.
    ///
    /// # Errors
    ///
    /// Returns a `ScrubbyError` if any download or unpacking operation fails.
    ///
    /// # Example
    ///
    /// ```
    /// let downloader = ScrubbyDownloader::new(outdir, indices);
    /// downloader.download_index();
    /// ```
    pub fn download_index(&self) -> Result<(), ScrubbyError> {

        if self.indices.is_empty() {
            log::warn!("No index names provided for download")
        }
        
        for index in &self.indices {
            for aligner in &self.aligners {
                let file_path = self.outdir.join(index.aligner_name(&aligner));
                log::info!("Downloading alignment index to file: {}", file_path.display());
                self.download(&index.aligner_name(aligner), &file_path)?;
                log::info!("Unpacking alignment index to directory: {}", self.outdir.display());
                self.unpack(&file_path, &self.outdir)?;
                log::info!("Removing download: {}", file_path.display());
                remove_file(&file_path)?;
            }
            for classifier in &self.classifiers {
                let file_path = self.outdir.join(index.classifier_name(&classifier));
                log::info!("Downloading classifier index to file: {}", file_path.display());
                self.download(&index.classifier_name(classifier), &file_path)?;
                log::info!("Unpacking classifier index to directory: {}", self.outdir.display());
                self.unpack(&file_path, &self.outdir)?;
                log::info!("Removing download: {}", file_path.display());
                remove_file(&file_path)?;
            }
        }

        Ok(())
    }
     /// Unpacks the downloaded file to the specified output directory.
    ///
    /// # Arguments
    ///
    /// * `download` - Path to the downloaded file.
    /// * `outdir` - Output directory for unpacked files.
    ///
    /// # Errors
    ///
    /// Returns a `ScrubbyError` if the unpacking operation fails.
    ///
    /// # Example
    ///
    /// ```
    /// let downloader = ScrubbyDownloader::new(outdir, indices);
    /// downloader.unpack(&download_path, &outdir);
    /// ```
    pub fn unpack(&self, download: &PathBuf, outdir: &PathBuf) -> Result<(), ScrubbyError> {
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
    /// Downloads a file from the specified URL to the given path.
    ///
    /// # Arguments
    ///
    /// * `file_name` - The name of the file to download.
    /// * `path` - The path where the file should be saved.
    ///
    /// # Errors
    ///
    /// Returns a `ScrubbyError` if the download operation fails.
    ///
    /// # Example
    ///
    /// ```
    /// let downloader = ScrubbyDownloader::new(outdir, indices);
    /// downloader.download("file_name.tar.xz", &path);
    /// ```
    pub fn download(&self, file_name: &str, path: &PathBuf) -> Result<(), ScrubbyError> {
        let url = format!("{}/{}", self.base_url, file_name);

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

/// Builder for creating an instance of `ScrubbyDownloader`.
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
    /// Creates a new instance of ScrubbyDownloaderBuilder.
    ///
    /// # Arguments
    ///
    /// * `outdir` - Output directory for downloaded files.
    /// * `indices` - A list of `ScrubbyIndex` to download.
    ///
    /// # Example
    ///
    /// ```
    /// let outdir = PathBuf::from("/path/to/output");
    /// let indices = vec![ScrubbyIndex::Chm13v2];
    /// let builder = ScrubbyDownloaderBuilder::new(outdir, indices);
    /// ```
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
    /// Sets the aligners for the builder.
    ///
    /// # Arguments
    ///
    /// * `aligner` - A list of `Aligner` instances.
    ///
    /// # Example
    ///
    /// ```
    /// let builder = ScrubbyDownloaderBuilder::new(outdir, indices).aligner(vec![Aligner::Bowtie2]);
    /// ```
    pub fn aligner<T: Into<Option<Vec<Aligner>>>>(mut self, aligner: T) -> Self {
        self.aligners = aligner.into();
        self
    }
    /// Sets the classifiers for the builder.
    ///
    /// # Arguments
    ///
    /// * `classifier` - A list of `Classifier` instances.
    ///
    /// # Example
    ///
    /// ```
    /// let builder = ScrubbyDownloaderBuilder::new(outdir, indices).classifier(vec![Classifier::new()]);
    /// ```
    pub fn classifier<T: Into<Option<Vec<Classifier>>>>(mut self, classifier: T) -> Self {
        self.classifiers = classifier.into();
        self
    }
    /// Sets the timeout duration for the builder.
    ///
    /// # Arguments
    ///
    /// * `timeout` - Timeout duration in seconds.
    ///
    /// # Example
    ///
    /// ```
    /// let builder = ScrubbyDownloaderBuilder::new(outdir, indices).timeout(60);
    /// ```
    pub fn timeout<T: Into<Option<u64>>>(mut self, timeout: T) -> Self {
        self.timeout = timeout.into();
        self
    }
    /// Sets the base URL for the builder.
    ///
    /// # Arguments
    ///
    /// * `base_url` - The base URL for downloading files.
    ///
    /// # Example
    ///
    /// ```
    /// let builder = ScrubbyDownloaderBuilder::new(outdir, indices).base_url("https://example.com");
    /// ```
    pub fn base_url<T: Into<Option<String>>>(mut self, base_url: T) -> Self {
        self.base_url = base_url.into();
        self
    }
    /// Sets the username for authentication.
    ///
    /// # Arguments
    ///
    /// * `username` - The username for basic authentication.
    ///
    /// # Example
    ///
    /// ```
    /// let builder = ScrubbyDownloaderBuilder::new(outdir, indices).username("user");
    /// ```
    pub fn username<T: Into<Option<String>>>(mut self, username: T) -> Self {
        self.username = username.into();
        self
    }
    /// Sets the password for authentication.
    ///
    /// # Arguments
    ///
    /// * `password` - The password for basic authentication.
    ///
    /// # Example
    ///
    /// ```
    /// let builder = ScrubbyDownloaderBuilder::new(outdir, indices).password("pass");
    /// ```
    pub fn password<T: Into<Option<String>>>(mut self, password: T) -> Self {
        self.password = password.into();
        self
    }
    /// Builds the `ScrubbyDownloader` instance.
    ///
    /// # Errors
    ///
    /// Returns a `ScrubbyError` if the downloader cannot be built.
    ///
    /// # Example
    ///
    /// ```
    /// let builder = ScrubbyDownloaderBuilder::new(outdir, indices);
    /// let downloader = builder.timeout(60).build();
    /// ```
    pub fn build(self) -> Result<ScrubbyDownloader, ScrubbyError> {
        
        if !self.outdir.exists() || !self.outdir.is_dir() {
            create_dir_all(&self.outdir)?;
        }
        
        let username = self.username
            .unwrap_or("u416706-sub1".to_string());
        let password = self.password
            .unwrap_or("G8tGWjHvdhUg4NGN".to_string());
        let base_url = self.base_url
            .unwrap_or(format!("https://{username}.your-storagebox.de/databases"));
        let aligners = self.aligners
            .unwrap_or(Vec::from([Aligner::Bowtie2]));
        let classifiers = self.classifiers
            .unwrap_or(Vec::new());
        let timeout = self.timeout
            .unwrap_or(30);

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
