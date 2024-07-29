use csv::WriterBuilder;
use env_logger::Builder;
use log::LevelFilter;
use needletail::{parse_fastx_file, FastxReader};
use niffler::get_writer;
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

/// Utility function to get a Needletail reader and Niffler compressed/uncompressed writer.
///
/// # Arguments
///
/// * `input` - A reference to the input file path.
/// * `output` - A reference to the output file path.
/// * `compression_level` - The desired compression level.
/// * `output_format` - Optional output format.
///
/// # Returns
///
/// * `Result<(Box<dyn FastxReader>, Box<dyn std::io::Write>), ScrubbyError>` - A tuple of reader and writer on success, otherwise an error.
///
/// # Example
///
/// ```
/// let (reader, writer) = get_niffler_fastx_reader_writer(&input_path, &output_path, niffler::compression::Level::Six, None).unwrap();
/// ```
pub fn get_niffler_fastx_reader_writer(
    input: &PathBuf,
    output: &PathBuf,
    compression_level: niffler::compression::Level,
    output_format: Option<niffler::compression::Format>,
) -> Result<(Box<dyn FastxReader>, Box<dyn std::io::Write>), ScrubbyError> {
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

pub struct DifferenceResult {
    pub input: u64,
    pub output: u64,
    pub difference: u64,
    pub read_ids: HashSet<String>
}


pub struct ReadDifference {
    pub result: DifferenceResult
}
impl ReadDifference {
    pub fn from(input_reads: &Vec<PathBuf>, output_reads: &Vec<PathBuf>) -> Result<Self, ScrubbyError> {
        let mut diff_ids = HashSet::new();

        let mut input_total = 0;
        let mut output_total = 0;
        let mut diff_total = 0;
        for (fq1, fq2) in input_reads.iter().zip(output_reads.iter()) {
            let mut reads2_ids = HashSet::new();

            let mut reader2 = parse_fastx_file(fq2)?;
            while let Some(record) = reader2.next() {
                let rec = record?;
                let read_id = get_id(&rec.id())?;
                reads2_ids.insert(read_id);
                output_total += 1;
            }
            let mut reader1 = parse_fastx_file(fq1)?;
            while let Some(record) = reader1.next() {
                let rec = record?;
                let read_id = get_id(&rec.id())?;
                if !reads2_ids.contains(&read_id) {
                    diff_ids.insert(read_id);
                    diff_total += 1
                }
                input_total += 1;
            }
        }
        Ok(Self {
            result: DifferenceResult { input: input_total, output: output_total, difference: diff_total, read_ids: diff_ids }
        })
    }
    pub fn write_reads(&self, output: &PathBuf, header: bool) -> Result<(), ScrubbyError> {
        
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

        for id in &self.result.read_ids {
            csv_writer.serialize(DifferenceRead { id: id.to_owned() })?;
        }

        csv_writer.flush()?;

        Ok(())
    }
}