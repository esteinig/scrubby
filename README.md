# scrubby <a href='https://github.com/esteinig'><img src='docs/scrubby.png' align="right" height="200" /></a>

[![build](https://github.com/esteinig/nanoq/actions/workflows/rust-ci.yaml/badge.svg?branch=master)](https://github.com/esteinig/scrubby/actions/workflows/rust-ci.yaml)
![](https://img.shields.io/badge/version-0.3.0-black.svg)

Host background depletion for metagenomic diagnostics with benchmarks and optimisation for clinical sequencing protocols and application scenarios.

## Overview

**`v1.0.0`**

- [Purpose](#purpose)
- [Install](#install)
- [Command-line interface](#command-line-interface)
- [Commands and options](#commands-and-options)
- [Rust library](#rust-library)
- [Dependencies](#dependencies)

## Purpose

## Install

Scrubby is available as binary release for MacOS and Linux; the default version requires several dependencies (aligners, classifiers) while the `mm2` release comes with an integrated and multi-threaded `minimap2-rs` implementation.

### Source

```
git clone https://github.com/esteinig/scrubby && cd scrubby
```

Compile default version with dependencies for aligners (+ samtools) or classifiers used:

```
cargo build --release
```

Compile built-in `minimap2-rs` version with `mm2` feature flag (experimental):

```
cargo build --release --features mm2
```

## Command-line interface

- Reads should be quality- and adapter-trimmed before applying `Scrubby`.
- Single or paired-end reads are supported with optional `gz` (`-i r1.fq r2.fq -o c1.fq.gz c2.fq.gz`). 
- Paired-end reads are always depleted/extracted as a pair (no unpaired read output).
- Default `minimap2` presets are `sr` for paired-end reads and `map-ont` for single reads.
- Multiple values can be specified consecutively or using multiple arguments (`-T Metazoa -T Bacteria`)


### Reference indices

List pre-built index names:

```
scrubby download --list
```

Download pre-built index by name for default aligner:

```
scrubby download --name chm13v2 --outdir . --aligner minimap2 bowtie2 --classifier kraken2
```

More options for aligners and classifier index download:

```
scrubby download --help
```

### Read depletion or extraction

Read depletion pipeline with `Bowtie2` aligner (default for paired-end reads):

```
scrubby reads -i R1.fq R2.fq -o C1.fq C2.fq -I chm13v2
```

Use built-in `minimap2-rs` if compiled with `mm2` feature (default for paired-end and long reads):

```
scrubby reads -i R1.fq R2.fq -o C1.fq C2.fq -I chm13v2.fa.gz
```

Long reads with non-default preset and `minimap2` aligner (default for logn reads):

```
scrubby reads -i R.fq.gz -o C.fq.gz -I chm13v2.fa.gz --preset lr-hq
```

Use classifier `Kraken2` or `Metabuli` instead of aligner:

```
scrubby reads -i R1.fq R2.fq -o C1.fq C2.fq -T Chordata -D 9606 -I chm13v2_k2/ -c kraken2
```

Use different aligner `strobealign` or `minimap2`:

```
scrubby reads -i R1.fq R2.fq -o C1.fq C2.fq -I chm13v2.fa.gz -a strobealign
```

With report output and depleted read identifiers:

```
scrubby reads -i R1.fq R2.fq -o C1.fq C2.fq -I chm13v2 -j report.json -r reads.tsv
```

Input and output compressed reads, increase threads and set working directory:

```
scrubby reads -i R1.fq.gz R2.fq.gz -o C1.fq.gz C2.fq.gz -I chm13v2 -w /tmp -t 16
```


### Read depletion or extraction from outputs

Classifier output cleaning (Kraken2, Metabuli):

```
scrubby classifier \
  --input R1.fq R2.fq \
  --output C1.fq C2.fq\
  --report kraken2.report \
  --reads kraken2.reads \
  --taxa Chordata \
  --taxa-direct 9606
```

Alignment output cleaning (.sam|.bam|.cram|.paf) or read identifier list (.txt). Alignment format is recognized from file extension or can be explicitly set with `--format`:

```
scrubby alignment  \
  --input R1.fq R2.fq \
  --output C1.fq C2.fq\
  --alignment alignment.paf \
  --min-len 50 \
  --min-cov 0.5 \
  --min-mapq 50 \
  --format paf
```

Add the `--extract` (`-e`) flag to any of the above tasks to reverse read depletion for read extraction:

```
scrubby reads --extract ...
```

Difference between input and output reads with optional counts and read identifier summaries:

```
scrubby diff -i R1.fq R2.fq -o C1.fq C2.fq -j counts.json -r reads.tsv
```


### Report output format

```json
{
  "version": "0.7.0",
  "date": "2024-07-30T06:50:15Z",
  "command": "scrubby reads -i smoke_R1.fastq.gz -i smoke_R2.fastq.gz -o test_R1.fq.gz -o test_R2.fq.gz --index /data/opt/scrubby_indices/chm13v2 --threads 16 --workdir /tmp/test --json test.json",
  "input": [
    "smoke_R1.fastq.gz",
    "smoke_R2.fastq.gz"
  ],
  "output": [
    "test_R1.fq.gz",
    "test_R2.fq.gz"
  ],
  "reads_in": 6678,
  "reads_out": 3346,
  "reads_removed": 3332,
  "reads_extracted": 0,
  "settings": {
    "aligner": "bowtie2",
    "classifier": null,
    "index": "/data/opt/scrubby_indices/chm13v2",
    "alignment": null,
    "reads": null,
    "report": null,
    "taxa": [],
    "taxa_direct": [],
    "classifier_args": null,
    "aligner_args": null,
    "min_len": 0,
    "min_cov": 0.0,
    "min_mapq": 0,
    "extract": false
  }
}
```

In this example, the `settings.aligner` is `null` if a `--classifier` is set.

## Commands and options

### Global options and commands

```bash
scrubby 0.7.0 
Eike Steinig (@esteinig)

Taxonomic read depletion for clinical metagenomic diagnostics

Usage: scrubby [OPTIONS] <COMMAND>

Commands:
  reads       Deplete or extract reads using aligners or classifiers
  classifier  Deplete or extract reads from classifier outputs (Kraken2/Metabuli)
  alignment   Deplete or extract reads from aligner output with additional filters (SAM/BAM/PAF/TXT)
  download    List available indices and download files for aligners and classfiers
  diff        Get read counts and identifiers of the difference between input and output read files
  help        Print this message or the help of the given subcommand(s)

Options:
  -l, --log-file <LOG_FILE>  Output logs to file instead of terminal
  -h, --help                 Print help (see more with '--help')
  -V, --version              Print version
```

### Pre-built reference downloads

```bash
List available indices and download files for aligners and classfiers

Usage: scrubby download [OPTIONS] --name [<NAME>...]

Options:
  -n, --name [<NAME>...]            Index name to download [possible values: chm13v2]
  -o, --outdir <OUTDIR>             Output directory for index download [default: .]
  -a, --aligner [<ALIGNER>...]      Download index for one or more aligners [possible values: bowtie2, minimap2, strobealign]
  -c, --classfier [<CLASSFIER>...]  Download index for one or more classifiers [possible values: kraken2, metabuli]
  -l, --list                        List available index names and exit
  -t, --timeout <TIMEOUT>           Download timeout in minutes - increase for large files and slow connections [default: 360]
  -h, --help                        Print help (see more with '--help')
```


### Read depletion or extraction

```bash
Deplete or extract reads using aligners or classifiers

Usage: scrubby reads [OPTIONS] --index <INDEX>

Options:
  -i, --input [<INPUT>...]              Input read files (optional .gz)
  -o, --output [<OUTPUT>...]            Output read files (optional .gz)
  -I, --index <INDEX>                   Reference index for aligner or classifier
  -e, --extract                         Read extraction instead of depletion
  -a, --aligner <ALIGNER>               Aligner to use, default is: Bowtie2 [possible values: bowtie2, minimap2, strobealign, minimap2-rs]
  -c, --classifier <CLASSIFIER>         Classifier to use [possible values: kraken2, metabuli]
  -T, --taxa [<TAXA>...]                Taxa and all sub-taxa to deplete using classifiers
  -D, --taxa-direct [<TAXA_DIRECT>...]  Taxa to deplete directly using classifiers
  -t, --threads <THREADS>               Number of threads to use for aligner and classifier [default: 4]
  -j, --json <JSON>                     Summary output file (.json)
  -w, --workdir <WORKDIR>               Optional working directory
  -r, --read-ids <READ_IDS>             Read identifier file (.tsv)
  -h, --help                            Print help (see more with '--help')
```

### Classifier outputs


```bash
Deplete or extract reads from classifier outputs (Kraken2, Metabuli)

Usage: scrubby classifier [OPTIONS] --report <REPORT> --reads <READS> --classifier <CLASSIFIER>

Options:
  -i, --input [<INPUT>...]              Input read files (optional .gz)
  -o, --output [<OUTPUT>...]            Output read files (optional .gz)
  -e, --extract                         Read extraction instead of depletion
  -k, --report <REPORT>                 Kraken-style report output from classifier
  -j, --reads <READS>                   Kraken-style read classification output
  -c, --classifier <CLASSIFIER>         Classifier output style [possible values: kraken2, metabuli]
  -T, --taxa [<TAXA>...]                Taxa and all sub-taxa to deplete using classifiers
  -D, --taxa-direct [<TAXA_DIRECT>...]  Taxa to deplete directly using classifiers
  -j, --json <JSON>                     Summary output file (.json)
  -w, --workdir <WORKDIR>               Optional working directory
  -r, --read-ids <READ_IDS>             Read identifier file (.tsv)
  -h, --help                            Print help (see more with '--help')
```


### Alignment outputs


```bash
Deplete or extract reads from aligner output with additional filters (SAM/BAM/PAF)

Usage: scrubby alignment [OPTIONS] --alignment <ALIGNMENT>

Options:
  -i, --input [<INPUT>...]     Input read files (can be compressed with .gz)
  -o, --output [<OUTPUT>...]   Output read files (can be compressed with .gz)
  -e, --extract                Read extraction instead of depletion
  -a, --alignment <ALIGNMENT>  Alignment file in SAM/BAM/PAF/TXT format
  -f, --format <FORMAT>        Explicit alignment format [possible values: sam, bam, cram, paf, txt]
  -l, --min-len <MIN_LEN>      Minimum query alignment length filter [default: 0]
  -c, --min-cov <MIN_COV>      Minimum query alignment coverage filter [default: 0]
  -q, --min-mapq <MIN_MAPQ>    Minimum mapping quality filter [default: 0]
  -j, --json <JSON>            Summary output file (.json)
  -w, --workdir <WORKDIR>      Optional working directory
  -r, --read-ids <READ_IDS>    Read identifier file (.tsv)
  -h, --help                   Print help (see more with '--help')
```

### Read difference

```bash
Get read counts and identifiers of the difference between input and output read files

Usage: scrubby diff [OPTIONS]

Options:
  -i, --input [<INPUT>...]    Input read files (.gz | .xz | .bz)
  -o, --output [<OUTPUT>...]  Output read files (.gz | .xz | .bz)
  -j, --json <JSON>           Summary output file (.json)
  -r, --read-ids <READ_IDS>   Read identifier file (.tsv)
  -h, --help                  Print help (see more with '--help')
```

## Rust library

You can use Scrubby with the builder structs from the prelude:

```rust
use scrubby::prelude::*;

let scrubby_mm2_ont = Scrubby::builder(
  "reads_in.fastq", 
  "reads_out.fastq"
)
  .json("report.json")
  .workdir("/tmp")
  .extract(false)
  .threads(16)
  .index("/path/to/reference.fasta")
  .aligner(Aligner::Minimap2)
  .preset(Preset::MapOnt)
  .build();

scrubby_mm2_ont.clean();

let scrubby_kraken2_metazoa = Scrubby::builder(
  "reads_in.fastq", 
  "reads_out.fastq"
)
  .json("report.json")
  .workdir("/tmp")
  .extract(false)
  .threads(16)
  .index("/path/to/kraken/index")
  .classifier(Classifier::Kraken2)
  .taxa(vec!["Metazoa"])
  .build();

scrubby_kraken2_metazoa.clean()

let scrubby_dl = ScrubbyDownloader::builder(
  "/path/to/download/directory", 
  vec![ScrubbyIndex::Chm13v2],
)
  .aligners(
    vec![Aligner::Minimap2, Aligner::Bowtie2]
  )
  .classifiers(
    vec![Classifier::Kraken2, Classifier::Metabuli]
  )
  .timeout(180)

scrubby_dl.list()
scrubby_dl.download_index();
```

## Dependencies

`Scrubby` wraps or implements the following libraries and tools:

* [`niffler`](https://github.com/luizirber/niffler)
* [`needletail`](https://github.com/onecodex/needletail)
* [`rust-htslib`](https://github.com/rust-bio/rust-htslib)
* [`minimap2`](https://github.com/lh3/minimap2)
* [`minimap2-rs`](https://github.com/jguhlin/minimap2-rs)
* [`strobealign`](https://github.com/ksahlin/strobealign)
* [`Bowtie2`](https://github.com/BenLangmead/bowtie2)
* [`Kraken2`](https://github.com/DerrickWood/kraken2)
* [`Metabuli`](https://github.com/steineggerlab/Metabuli)
