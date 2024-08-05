import typer 

from typing import List

from tensorflow.config import list_physical_devices
from tensorflow.test import is_built_with_cuda

from .model import train_nn, predict_nn

app = typer.Typer()


@app.command()
def train(
    fastq_files: List[str] = typer.Argument(..., help="Paths to the input FASTQ files."),
    model_weights: str = typer.Argument(..., help="Path to save the model weights."),
    alignment_data: str = typer.Option(None, help="Path to alignment data."),
    epochs: int = typer.Option(10, help="Number of epochs for training."),
    batch_size: int = typer.Option(32, help="Batch size for training."),
    multi_gpu: bool =  typer.Option(False, help="Multiple GPUs for training."),
):
    """Train the neural network."""
    train_nn(fastq_files, model_weights, alignment_data, epochs, batch_size)

@app.command()
def predict(
    model_weights: str = typer.Argument(..., help="Path to the model weights."),
    fastq_files: List[str] = typer.Argument(..., help="Paths to the input FASTQ files."),
    alignment_data: str = typer.Option(None, help="Path to alignment data.")
):
    """Predict classes for the input sequences."""
    predict_nn(model_weights, fastq_files, alignment_data)

@app.command()
def check_gpu(): 
    
    """ Check if Tensorflow detects GPU """
    # Check if TensorFlow is built with GPU support
    print("Is TensorFlow built with GPU support:", is_built_with_cuda())

    # List available GPUs
    gpus = list_physical_devices('GPU')
    print("Available GPUs:", gpus)

    # Verify that TensorFlow is using the GPU
    if gpus:
        print("TensorFlow is using the GPU")
    else:
        print("TensorFlow is not using the GPU")