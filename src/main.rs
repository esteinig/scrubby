use clap::Parser;
use scrubby::identity::{train_nn, predict_nn, check_gpu_connectivity};
use scrubby::utils::init_logger;
use scrubby::terminal::{App, Commands};


fn main() -> anyhow::Result<()> {
    
    let cli = App::parse();

    init_logger(cli.log_file);

    match cli.command {
        Commands::Reads(args) => {
            args.validate_and_build()?.clean()?;
        },
        Commands::Classifier(args) => {
            args.validate_and_build()?.clean()?;
        },
        Commands::Alignment(args) => {
            args.validate_and_build()?.clean()?;
        },
        Commands::Download(args) => {
            let dl = args.clone().validate_and_build()?;

            if args.list { dl.list(); } else { dl.download_index()?; }
        },
        Commands::Diff(args) => {
            args.validate_and_build()?.compute()?;
        },
        Commands::Nn(args) => {
            if args.train { 
                train_nn(args.device, args.fastq, args.model_weights, args.alignment, args.epochs as i64, args.batch_size, 10000)?;
            } else if args.check {
                if check_gpu_connectivity() {
                    log::info!("Successfully connected to the GPU.");
                } else {
                    log::info!("Failed to connect to the GPU.");
                }
            } else {
                predict_nn(args.device, args.model_weights, args.fastq, args.alignment)?;
            }
        },
    }

    Ok(())
}