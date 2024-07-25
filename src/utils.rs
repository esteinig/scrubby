use env_logger::Builder;
use log::LevelFilter;
use std::{fs::OpenOptions, io::Write, path::PathBuf};
use termcolor::{Color, ColorChoice, ColorSpec, StandardStream, WriteColor};

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
            writeln!(
                &mut stderr,
                "{}",
                record.level()
            )
            .unwrap();

            stderr.set_color(&default_style).unwrap();
            write!(&mut stderr, "] - {}", record.args()).unwrap();

            stderr.reset().unwrap();
        } else {
            
            stdout.set_color(&default_style).unwrap();
            write!(&mut stdout, "{} [", timestamp).unwrap();

            stdout.set_color(&style).unwrap();
            writeln!(
                &mut stdout,
                "{}",
                record.level()
            )
            .unwrap();

            stdout.set_color(&default_style).unwrap();
            write!(&mut stdout, "] - {}", record.args()).unwrap();

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
