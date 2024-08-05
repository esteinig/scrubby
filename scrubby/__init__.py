import typer
from pathlib import Path

app = typer.Typer(add_completion=False)


@app.command()
def train(
    benchmark_dir: Path = typer.Option(
        ..., help="Benchmark workflow output directory from Cipher"
    ),
    output_pdf: Path = typer.Option(
        ..., help="Output file path (PDF)"
    ),
):
    """ Process host depletion evaluation benchmark """
