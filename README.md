# scrubby <a href='https://github.com/esteinig'><img src='docs/scrubby.png' align="right" height="200" /></a>

[![build](https://github.com/esteinig/nanoq/actions/workflows/rust-ci.yaml/badge.svg?branch=master)](https://github.com/esteinig/scrubby/actions/workflows/rust-ci.yaml)
![](https://img.shields.io/badge/version-1.0.0-black.svg)

Host background depletion for metagenomic diagnostics with benchmarks and optimisation for clinical sequencing protocols and application scenarios.

## Overview

- [Purpose](#purpose)
- [Install](#install)
- [Command-line interface](#command-line-interface)
- [Commands and options](#commands-and-options)
- [Rust library](#rust-library)
- [Dependencies](#dependencies)

## Purpose

## Install

Scrubby is available as binary release for Linux and macOS (x86_64); the default version requires several dependencies (aligners + `samtools` or classifiers used for read depletion or extraction); the `mm2` release comes with a multi-threaded `minimap2-rs` implementation and does not require additional dependencies.

### Source

```
git clone https://github.com/esteinig/scrubby && cd scrubby
```

Compile default version - requires aligners and classifiers (+ `samtools`):

```
cargo build --release
```

Compile built-in `minimap2-rs` version with `mm2` feature flag (experimental):

```
cargo build --release --features mm2
```

## Command-line interface

- Reads should be quality- and adapter-trimmed before applying `Scrubby`.
- Single or paired-end reads are supported with optional `gz` input/output compression. 
- Paired-end reads are always depleted/extracted as a pair (no unpaired read output).
- Default `minimap2` presets are `sr` for paired-end reads and `map-ont` for single reads.
- Multiple values can be specified consecutively or using multiple arguments (`-T Metazoa -T Bacteria`)


### Reference indices

List pre-built index names:

```shell
scrubby download --list
```

Download pre-built index by name for default aligner:

```shell
scrubby download --name chm13v2 --outdir . --aligner minimap2 bowtie2 --classifier kraken2
```

More options for aligners and classifier index download:

```shell
scrubby download --help
```

### Read depletion or extraction

Read depletion pipeline with `Bowtie2` aligner (default for paired-end reads):

```shell
scrubby reads -i r1.fq r2.fq -o c1.fq c2.fq -I chm13v2
```

Use built-in `minimap2-rs` if compiled with `mm2` feature (default for paired-end and long reads):

```shell
scrubby reads -i r1.fq r2.fq -o c1.fq c2.fq -I chm13v2.fa.gz
```

Long reads with non-default preset and `minimap2` aligner (default for long reads):

```shell
scrubby reads -i r.fq -o c.fq -I chm13v2.fa.gz --preset lr-hq
```

Single-end short reads requires explicit aligner and preset for `minimap2`:

```shell
scrubby reads -i r1.fq -o c1.fq -I chm13v2.fa.gz --aligner minimap2 --preset sr
```

Use classifier `Kraken2` or `Metabuli` instead of aligner:

```shell
scrubby reads -i r1.fq r2.fq -o c1.fq c2.fq -T Chordata -D 9606 -I chm13v2_k2/ -c kraken2
```

Use different aligner `strobealign` or `minimap2`:

```shell
scrubby reads -i r1.fq r2.fq -o c1.fq c2.fq -I chm13v2.fa.gz -a strobealign
```

With report output and depleted read identifiers:

```shell
scrubby reads -i r1.fq r2.fq -o c1.fq c2.fq -I chm13v2 -j report.json -r reads.tsv
```

Input and output compressed reads, increase threads and set working directory:

```shell
scrubby reads -i r1.fq.gz r2.fq.gz -o c1.fq.gz c2.fq.gz -I chm13v2 -w /tmp -t 16
```

### Read depletion or extraction directly from outputs

Classifier output cleaning for Kraken-style reports and read classification outputs (Kraken2, Metabuli):

```shell
scrubby classifier \
  --input r1.fq r2.fq \
  --output c1.fq c2.fq\
  --report kraken2.report \
  --reads kraken2.reads \
  --taxa Chordata \
  --taxa-direct 9606
```

Alignment output cleaning (.sam|.bam|.cram|.paf) or read identifier list (.txt). Alignment format is recognized from file extension or can be explicitly set with `--format`. Alignment can be '-' for reading from `stdin` with explicit format argument. PAF and TXT formats can be compressed (.gz|.xz|.bz) unless reading
from `stdin`.

```shell
scrubby alignment  \
  --input r1.fq r2.fq \
  --output c1.fq c2.fq\
  --alignment alignment.paf \
  --min-len 50 \
  --min-cov 0.5 \
  --min-mapq 50

minimap2 -x map-ont ref.fa r.fq | scrubby alignment -a - -f paf -i r.fq -o c.fq
```

Add the `--extract` (`-e`) flag to any of the above tasks to reverse read depletion for read extraction:

```shell
scrubby reads --extract ...
```

Difference between input and output reads with optional counts and read identifier summaries:

```shell
scrubby diff -i r1.fq r2.fq -o c1.fq c2.fq -j counts.json -r reads.tsv
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

```shell
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

```shell
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

```shell
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


```shell
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


```shell
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

```shell
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


// Example running Minimap2 on long reads

let scrubby_mm2_ont = Scrubby::builder(
  "/path/to/reads_in.fastq", 
  "/path/to/reads_out.fastq"
)
  .json("/path/to/report.json")
  .extract(false)
  .threads(16)
  .index("/path/to/reference.fasta")
  .aligner(Aligner::Minimap2)
  .preset(Preset::MapOnt)
  .build();

scrubby_mm2_ont.clean();


// Example running Minimap2 on paired-end reads

let scrubby_mm2_sr = Scrubby::builder(
  vec!["/path/to/reads_in_R1.fastq", "/path/to/reads_in_R2.fastq"] 
  vec!["/path/to/reads_out_R1.fastq", "/path/to/reads_out_R2.fastq"]
)
  .json("/path/to/report.json")
  .extract(false)
  .threads(16)
  .index("/path/to/reference.fasta")
  .aligner(Aligner::Minimap2)
  .preset(Preset::Sr)
  .build();

scrubby_mm2_sr.clean();


// Example running Kraken2, depleting Metazoa 

let scrubby_kraken2_metazoa = Scrubby::builder(
  "/path/to/reads_in.fastq", 
  "/path/to/reads_out.fastq"
)
  .json("/path/to/report.json")
  .extract(false)
  .threads(16)
  .index("/path/to/kraken/index")
  .classifier(Classifier::Kraken2)
  .taxa(vec!["Metazoa"])
  .build();

scrubby_kraken2_metazoa.clean();


// Example from Kraken2 outputs, depleting Metazoa 

let scrubby_kraken2_output_metazoa = Scrubby::builder(
  "/path/to/reads_in.fastq", 
  "/path/to/reads_out.fastq"
)
  .json("/path/to/report.json")
  .extract(false)
  .classifier(Classifier::Kraken2)
  .report("/path/to/kraken/report")
  .reads("/path/to/kraken/read/classifications")
  .taxa(vec!["Metazoa"])
  .build();

scrubby_kraken2_output_metazoa.clean();


// Example from alignment output file with filters

let scrubby_paf_output_filters = Scrubby::builder(
  "/path/to/reads_in.fastq", 
  "/path/to/reads_out.fastq"
)
  .json("/path/to/report.json")
  .extract(false)
  .alignment("/path/to/aln.paf")
  .min_query_length(50)
  .min_query_coverage(0.5)
  .min_mapq(50)
  .build();

scrubby_paf_output_filters.clean();

// Downloader example

let scrubby_dl = ScrubbyDownloader::builder(
  "/path/to/download/directory", 
  vec![ScrubbyIndex::Chm13v2]
)
  .timeout(180)
  .aligner(vec![
    Aligner::Minimap2, Aligner::Bowtie2
  ])
  .classifier(vec![
    Classifier::Kraken2, Classifier::Metabuli
  ])
  .build();

scrubby_dl.list();
scrubby_dl.download_index();
```

## Dependencies

Rust libraries:

* [`niffler`](https://github.com/luizirber/niffler)
* [`needletail`](https://github.com/onecodex/needletail)
* [`rust-htslib`](https://github.com/rust-bio/rust-htslib)
* [`minimap2-rs`](https://github.com/jguhlin/minimap2-rs)

Aligners and classifiers:

* [`samtools`](https://github.com/samtools/samtools)
* [`minimap2`](https://github.com/lh3/minimap2)
* [`strobealign`](https://github.com/ksahlin/strobealign)
* [`Bowtie2`](https://github.com/BenLangmead/bowtie2)
* [`Kraken2`](https://github.com/DerrickWood/kraken2)
* [`Metabuli`](https://github.com/steineggerlab/Metabuli)
