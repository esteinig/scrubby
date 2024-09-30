use csv::WriterBuilder;
use env_logger::Builder;
use log::LevelFilter;
use needletail::{parse_fastx_file, FastxReader};
use niffler::{get_reader, get_writer};
use serde::{Deserialize, Serialize};
use std::{collections::HashSet, ffi::OsStr, fs::{File, OpenOptions}, io::{BufWriter, Write}, path::{Path, PathBuf}};
use termcolor::{Color, ColorChoice, ColorSpec, StandardStream, WriteColor};

use crate::error::ScrubbyError;


/// Extension trait for inferring compression format from file extension.
pub trait CompressionExt {
    fn from_path<S: AsRef<OsStr> + ?Sized>(p: &S) -> Self;
}

impl CompressionExt for niffler::compression::Format {
    /// Attempts to infer the compression type from the file extension.
    /// If the extension is not known, `Uncompressed` is returned.
    ///
    /// # Example
    ///
    /// ```
    /// use scrubby::CompressionExt;
    /// let format = niffler::compression::Format::from_path("file.gz");
    /// ```
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

/// Utility function to get a Needletail and Niffler compressed/uncompressed writer.
///
/// # Arguments
///
/// * `output` - A reference to the output file path.
/// * `compression_level` - The desired compression level.
/// * `output_format` - Optional output format.
///
/// # Returns
///
/// * `Result<(Box<dyn FastxReader>, Box<dyn std::io::Write>), ScrubbyError>` - A writer on success, otherwise an error.
///
/// # Example
///
/// ```
/// let (reader, writer) = get_niffler_fastx_reader_writer(&output_path, niffler::compression::Level::Six, None).unwrap();
/// ```
pub fn get_fastx_writer(
    output: &PathBuf,
    compression_level: niffler::compression::Level,
    output_format: Option<niffler::compression::Format>,
) -> Result<Box<dyn std::io::Write>, ScrubbyError> {
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

    Ok(writer)
}

/// Utility function to extract the ID from a FASTQ record header.
///
/// # Arguments
///
/// * `id` - A byte slice containing the FASTQ record header.
///
/// # Returns
///
/// * `Result<String, ScrubbyError>` - The extracted ID as a string on success, otherwise an error.
///
/// # Example
///
/// ```
/// let id = get_id(b"@read1 description").unwrap();
/// ```
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


pub fn init_logger(log_file: Option<PathBuf>) {
    // Create a buffer for writing log messages with colors
    let mut builder = Builder::new();

    // Define the log format
    builder.format(|buf, record| {
        let mut stdout = StandardStream::stdout(ColorChoice::Always);
        let mut stderr = StandardStream::stderr(ColorChoice::Always);

        let mut style = ColorSpec::new();
        match record.level() {
            log::Level::Trace => style.set_fg(Some(Color::White)).set_bold(true),
            log::Level::Debug => style.set_fg(Some(Color::Rgb(255, 195, 0))).set_bold(true),
            log::Level::Info => style.set_fg(Some(Color::Green)).set_bold(true),
            log::Level::Warn => style.set_fg(Some(Color::Rgb(255, 102, 0))).set_bold(true),
            log::Level::Error => style.set_fg(Some(Color::Red)).set_bold(true),
        };

        let mut default_style = ColorSpec::new();
        default_style.set_fg(Some(Color::White));

        let timestamp = buf.timestamp();

        if record.level() == log::Level::Error || record.level() == log::Level::Warn {
            
            stderr.set_color(&default_style).unwrap();
            write!(&mut stderr, "{} [", timestamp).unwrap();

            stderr.set_color(&style).unwrap();
            write!(&mut stderr,"{}", record.level()).unwrap();

            stderr.set_color(&default_style).unwrap();
            writeln!(&mut stderr, "] - {}", record.args()).unwrap();

            stderr.reset().unwrap();
        } else {
            
            stdout.set_color(&default_style).unwrap();
            write!(&mut stdout, "{} [", timestamp).unwrap();

            stdout.set_color(&style).unwrap();
            write!(&mut stdout, "{}", record.level()).unwrap();

            stdout.set_color(&default_style).unwrap();
            writeln!(&mut stdout, "] - {}", record.args()).unwrap();

            stdout.reset().unwrap();
        }

        Ok(())
    });

    // Set the default filter level to Info
    builder.filter(None, LevelFilter::Info);

    if let Some(log_path) = log_file {
        let file = OpenOptions::new()
            .append(true)
            .create(true)
            .open(log_path)
            .unwrap();

        builder.target(env_logger::Target::Pipe(Box::new(file)));
    }

    builder.init();
}


#[derive(Serialize, Deserialize)]
pub struct DifferenceRead {
    id: String
}

#[derive(Serialize, Deserialize)]
pub struct Difference {
    pub reads_in: u64,
    pub reads_out: u64,
    pub difference: u64,
    #[serde(skip_serializing)]
    pub read_ids: HashSet<String>
}
impl Difference {
    pub fn new(reads_in: u64, reads_out: u64, difference: u64) -> Self {
        Self { reads_in, reads_out, difference, read_ids: HashSet::new() }
    }
    pub fn to_json(&self, output: &PathBuf) -> Result<(), ScrubbyError> {
        let mut file = std::fs::File::create(output)?;
        let json_string = serde_json::to_string_pretty(self)?;
        file.write_all(json_string.as_bytes())?;
        Ok(())
    }
    pub fn write_read_ids(&self, output: &PathBuf, header: bool) -> Result<(), ScrubbyError> {
        
        let buf_writer = BufWriter::new(File::create(&output)?);
        let writer = get_writer(
            Box::new(buf_writer), 
            niffler::Format::from_path(output), 
            niffler::compression::Level::Nine
        )?;

        let mut csv_writer = WriterBuilder::new()
            .delimiter(b'\t')
            .has_headers(header)
            .from_writer(writer);

        for id in &self.read_ids {
            csv_writer.serialize(DifferenceRead { id: id.to_owned() })?;
        }

        csv_writer.flush()?;

        Ok(())
    }
}


pub struct ReadDifference {
    pub input_reads: Vec<PathBuf>,
    pub output_reads: Vec<PathBuf>,
    pub json: Option<PathBuf>,
    pub read_ids: Option<PathBuf>
}
impl ReadDifference {
    pub fn new(input_reads: &Vec<PathBuf>, output_reads: &Vec<PathBuf>, json: Option<PathBuf>, read_ids: Option<PathBuf>) -> Self {
        Self {
            input_reads: input_reads.to_owned(),
            output_reads: output_reads.to_owned(),
            json,
            read_ids
        }
    }
    pub fn compute(&self) -> Result<Difference, ScrubbyError> {
        let diff = self.get_difference()?;
            
        if let Some(json) = &self.json {
            diff.to_json(json)?;
        }
        if let Some(read_ids) = &self.read_ids {
            diff.write_read_ids(read_ids, true)?;
        }

        Ok(diff)
    }
    pub fn get_difference(&self) -> Result<Difference, ScrubbyError> {
        let mut diff_ids = HashSet::new();

        let mut input_total = 0;
        let mut output_total = 0;
        let mut diff_total = 0;
        for (fq1, fq2) in self.input_reads.iter().zip(self.output_reads.iter()) {
            let mut reads2_ids = HashSet::new();

            let reader2 = parse_fastx_file_with_check(fq2)?;
            if let Some(mut reader2) = reader2 {
                while let Some(record) = reader2.next() {
                    let rec = record?;
                    let read_id = get_id(&rec.id())?;
                    reads2_ids.insert(read_id);
                    output_total += 1;
                }
            }
            
            let reader1 = parse_fastx_file_with_check(fq1)?;
            if let Some(mut reader1) = reader1 {
                while let Some(record) = reader1.next() {
                    let rec = record?;
                    let read_id = get_id(&rec.id())?;
                    if !reads2_ids.contains(&read_id) {
                        diff_ids.insert(read_id);
                        diff_total += 1
                    }
                    input_total += 1;
                }
            } else {
                log::warn!("Input file is empty: {}", fq1.display())
            }
        }
        Ok(Difference { reads_in: input_total, reads_out: output_total, difference: diff_total, read_ids: diff_ids })
    }
}


pub struct ReadDifferenceBuilder {
    input_reads: Vec<PathBuf>,
    output_reads: Vec<PathBuf>,
    read_ids: Option<PathBuf>,
    json: Option<PathBuf>
}

impl ReadDifferenceBuilder {
    pub fn new(input_reads: &Vec<PathBuf>, output_reads: &Vec<PathBuf>) -> Self {
        ReadDifferenceBuilder {
            input_reads: input_reads.clone(),
            output_reads: output_reads.clone(),
            read_ids: None,
            json: None
        }
    }
    /// Sets the `read_ids` field.
    ///
    /// # Example
    ///
    /// ```
    /// use scrubby::ReadDifferenceBuilder;
    /// use std::path::PathBuf;
    ///
    /// let builder = ReadDifferenceBuilder::new(...).read_ids(PathBuf::from("depleted_reads.fastq"));
    /// ```
    pub fn read_ids<T: Into<Option<PathBuf>>>(mut self, read_ids: T) -> Self {
        self.read_ids = read_ids.into();
        self
    }
    /// Sets the `json` field.
    ///
    /// # Example
    ///
    /// ```
    /// use scrubby::ReadDifferenceBuilder;
    /// use std::path::PathBuf;
    ///
    /// let builder = ReadDifferenceBuilder::new(...).json(PathBuf::from("report.json"));
    /// ```
    pub fn json<T: Into<Option<PathBuf>>>(mut self, json: T) -> Self {
        self.json = json.into();
        self
    }

    pub fn build(self) -> Result<ReadDifference, ScrubbyError> {
        
        // Check if input and output vectors are not empty
        if self.input_reads.is_empty() || self.output_reads.is_empty() {
            return Err(ScrubbyError::EmptyInputOutput);
        }
        // Check if input and output vectors have the same length
        if self.input_reads.len() != self.output_reads.len() {
            return Err(ScrubbyError::MismatchedInputOutputLength);
        }
        // Check if input and output vectors length is limited to one or two
        if self.input_reads.len() > 2 || self.output_reads.len() > 2 {
            return Err(ScrubbyError::InputOutputLengthExceeded);
        }
        // Check if each input file exists and is a file
        for input_file in &self.input_reads {
            if !input_file.exists() || !input_file.is_file() {
                return Err(ScrubbyError::MissingInputReadFile(input_file.clone()));
            }
        }
        
        Ok(ReadDifference::new(&self.input_reads, &self.output_reads, self.json, self.read_ids))
    }
}

pub fn is_file_empty<P: AsRef<Path>>(path: P) -> Result<bool, ScrubbyError> {
    let file = File::open(&path)?;
    
    // Use niffler to get a reader for the (possibly compressed) file
    let (mut reader, _format) = match get_reader(Box::new(file)) {
        Ok(reader_format) => reader_format,
        Err(niffler::Error::FileTooShort) => return Ok(true),
        Err(e) => return Err(ScrubbyError::NifflerError(e)),
    };
    // Try to read the first byte
    let mut buffer = [0; 1];
    match reader.read(&mut buffer) {
        Ok(0) => Ok(true),
        Ok(_) => Ok(false), // Successfully read a byte, file is not empty
        Err(e) => Err(ScrubbyError::IoError(e))
    }
}

pub fn parse_fastx_file_with_check<P: AsRef<Path>>(path: P) -> Result<Option<Box<dyn FastxReader>>, ScrubbyError> {
    if is_file_empty(&path)? {
        Ok(None)
    } else {
        Ok(Some(parse_fastx_file(&path)?))
    }
}



pub trait IntoVecPathBuf {
    fn into_vec_path_buf(self) -> Vec<PathBuf>;
}

impl IntoVecPathBuf for &str {
    fn into_vec_path_buf(self) -> Vec<PathBuf> {
        vec![PathBuf::from(self)]
    }
}

impl IntoVecPathBuf for String {
    fn into_vec_path_buf(self) -> Vec<PathBuf> {
        vec![PathBuf::from(self)]
    }
}

impl IntoVecPathBuf for Vec<&str> {
    fn into_vec_path_buf(self) -> Vec<PathBuf> {
        self.into_iter().map(PathBuf::from).collect()
    }
}

impl IntoVecPathBuf for Vec<String> {
    fn into_vec_path_buf(self) -> Vec<PathBuf> {
        self.into_iter().map(PathBuf::from).collect()
    }
}

impl IntoVecPathBuf for Vec<PathBuf> {
    fn into_vec_path_buf(self) -> Vec<PathBuf> {
        self
    }
}