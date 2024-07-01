# scrubby <a href='https://github.com/esteinig'><img src='docs/scrubby.png' align="right" height="200" /></a>

[![build](https://github.com/esteinig/nanoq/actions/workflows/rust-ci.yaml/badge.svg?branch=master)](https://github.com/esteinig/scrubby/actions/workflows/rust-ci.yaml)
![](https://img.shields.io/badge/version-0.3.0-black.svg)

A (t)rusty read scrubber to deplete/extract background taxa using k-mer classifications or alignments. 

## Overview

**`v0.3.0`**

- [Purpose](#purpose)
- [Install](#install)
- [Usage](#usage)
  - [General options](#general-options)
  - [Read scrubbing](#read-scrubbing)
    - [Read scrubbing pipeline](#read-scrubbing-pipeline)
    - [Scrub reads with Kraken2](#kraken2-scrubbing)
    - [Scrub reads with alignments](#alignment-scrubbing)
    - [Summary report output](#summary-output)
- [Command-line arguments](#command-line-arguments)
  - [Read scrubbing pipeline](#read-scrubbing-pipeline)
  - [Database scrubbing pipeline](#database-scrubbing-pipeline)
  - [Kraken2 depletion/extraction](#kraken-scrubber)
  - [Alignment depletion/extraction](#alignment-scrubber)
- [Considerations](#considerations)
  - [Taxonomic database errors](#taxonomic-database-errors)
- [Roadmap](#roadmap)
- [Dependencies](#dependencies)

## Purpose

`Scrubby` can deplete/extract reads classified at taxonomic sub-ranks and perform sequential depletion/extraction from multiple databases or alignments. 

As an example, you can specify a primary (fast) k-mer depletion of all reads classified as Eukaryota (including sub-ranks like Holozoa) with `Kraken2`, then follow up with `minimap2` alignment against [`CHM13v2`](https://github.com/marbl/CHM13) to further deplete pesky human reads.

`Scrubby` is the mirror counterpart of the excellent [`ReadItAndKeep`](https://github.com/GlobalPathogenAnalysisService/read-it-and-keep) and meant to facilitate safe and thorough background depletion of host and other taxa. You can also use `Scrubby` for target retention by specifying the target reference/taxon and `--extract` flag to retain all reads aligned/classified by one or multiple methods.

This is a preliminary release, use at your own peril :skull:

## Install

Development version (v0.3.0)

```
git clone https://github.com/esteinig/scrubby
cd scrubby && cargo build --release
./target/release/scrubby --help
```

BioConda:

```
conda install -c conda-forge -c bioconda scrubby
```

Cargo:

```
cargo install scrubby
```

## Usage

Scrubbing pipeline:

```
scrubby scrub-reads --help
```

Scrubbing `Kraken2`:

```
scrubby scrub-kraken --help
```

Scrubbing alignment:

```
scrubby scrub-alignment --help
```

Add the `--extract` flag to any of the above tasks to enable read extraction:

```
scrubby scrub-reads --extract ...
```

### General options

Reads:

- Reads should be quality- and adapter-trimmed before applying `Scrubby`.
- Single or paired-end reads are supported (`--input r1.fq r2.fq --output c1.fq c2.fq`). 
- Paired-end reads are always depleted/extracted as a pair (no unpaired read output).
- Compression formats are recognized from extensions of `--input/--output` (`gz|bz|bz2|xz`).

Filters:

- Taxa for `Kraken2`/`Metabuli` can be `taxids` or `names` as listed in the report file (case sensitive).
- Alignment filters as in `ReadItAndKeep` can be specified (`--min-len`, `--min-cov`, `--min-mapq`). 
- Read depletion/extraction summaries can be written to file (`--json file.json`) or stdout (`--json -`). 
- Arguments for which multiple values can be supplied e.g. inputs/outputs (`-i/-o`), databases/references (`-k/-m/-b/-s`) or taxa (`-t/-d`) can be specified either consecutively (e.g. `-k Metazoa Bacteria`) or using multiple arguments (e.g. `-k Metazoa -k Bacteria`)

### Read scrubbing

#### Read scrubbing pipeline


`Scrubby` primarily depletes/extracts using sequential k-mer and alignment methods. This will call `Kraken2` and aligners (`minimap2`, `bowtie2`, `strobealign`) under the hood,
creating intermediary files in the `-W/--wordir` which can be retained (`-K/--keep`). By default the working directory is created with a time-stamp (`Scrubby_{YYYYMMDDTHHMMSS}`).


```
scrubby scrub-reads \
  --input R1.fq.gz R2.fq.gz \
  --output S1.fq.gz S2.fq.gz \
  --kraken-db minikraken/ \
  --kraken-taxa Metazoa \
  --minimap2-index chm13v2.fasta \
  --min-len 50
```


When using `Kraken2` reads classified within a particular taxonomic rank are depleted/extracted **except those above Domain**. Minor rank designations are considered to be sub-ranks (D1, D2, ...). If you only want to deplete/extract a specific taxon or want to deplete/extract a rank above Domain use the `--kraken-taxa-direct` argument (for example "Unclassified" and "Cellular Organisms") :

```
scrubby scrub-reads \
  --input R1.fq.gz R2.fq.gz \
  --output S1.fq.gz S2.fq.gz \
  --kraken-db minikraken/ \
  --kraken-taxa Eukaryota \
  --kraken-taxa-direct Unclassified 131567
```

You can skip methods by leaving out `--kraken-db` or `--minimap2-index` arguments, for example:

```
scrubby scrub-reads \
  --input R1.fq.gz R2.fq.gz \
  --output S1.fq.gz S2.fq.gz \
  --kraken-db minikraken/ \
  --kraken-taxa Eukaryota
```

Order of databases/references used for depletion is the same as the order of the command-line arguments. Order of methods of depletion is 
currently always: `Kraken2` --> `minimap2` --> `bowtie2` --> `strobealign` (unless a method is not specified).

#### Kraken2 scrubbing

You can deplete/extract a set of reads using pre-computed outputs from `Kraken2`, which requires the `--kraken-report` and `--kraken-reads` files:

```
scrubby scrub-reads \
  --input R1.fq.gz R2.fq.gz \
  --output S1.fq.gz S2.fq.gz \
  --kraken-report test.report \
  --kraken-reads test.reads \
  --kraken-taxa Metazoa
```

#### Alignment scrubbing

You can deplete/extract a set of reads using pre-computed alignments in `PAF|SAM|BAM|CRAM` formats. These are recognized from the alignment extension 
or can be specified explicitly with `--alignment-format`:

```
scrubby scrub-reads \
  --input R1.fq.gz R2.fq.gz \
  --output S1.fq.gz S2.fq.gz \
  --kraken-report test.report \
  --kraken-reads test.reads \
  --kraken-taxa Metazoa
```

#### Summary output

The schema contains a `pipeline` array for each database or reference provided in the order in which reads were depleted/extracted in the read scrubbing pipeline. Tool values are lowercase tool names, one of: `kraken2`, `minimap2`, `strobealign`. 

Note that when individually depleting/extracting `Kraken2` (`scrubby scrub-kraken`) or alignments (`scrubby scrub-alignments`):

  - the `pipeline` array contains a single entry which always has an index of `0`
  - the `path` value is the path to the classified reads file for `Kraken2`
  - the `tool` value in the pipeline array entry is `null`

```json
{
  "version": "0.3.0",
  "schema_version": "0.3.0",
  "settings": {
    "kraken_taxa": [
      "Eukaryota",
      "Bacteria"
    ],
    "kraken_taxa_direct": [
      "Unclassified"
    ],
    "min_len": 30,
    "min_cov": 0.0,
    "min_mapq": 0,
    "extract": false
  },
  "summary": {
    "total": 1000000,
    "depleted": 66784,
    "extracted": 0
  },
  "pipeline": [
    {
      "index": 0,
      "tool": "kraken2",
      "name": "SILVA_138_rRNA",
      "path": "/path/to/SILVA_138_rRNA",
      "total": 1000000,
      "depleted": 66784,
      "extracted": 0,
      "files": [
        {
          "total": 500000,
          "depleted": 33392,
          "extracted": 0,
          "input_file": "/path/to/test_r1.fq",
          "output_file": "/path/to/tmp/workdir/0-rrna_1.fq"
        },
        {
          "total": 500000,
          "depleted": 33392,
          "extracted": 0,
          "input_file": "/path/to/test_r2.fq",
          "output_file": "/path/to/tmp/workdir/0-rrna_2.fq"
        }
      ]
    }
  ]
}
```

## Command-line arguments

### Scrubbing pipeline

```shell
scrubby-scrub-reads 0.3.0
Clean sequence reads by removing background taxa (Kraken2) or aligning reads (Minimap2)

USAGE:
    scrubby scrub-reads [FLAGS] [OPTIONS] --input <input>... --kraken-db <kraken-db>... --output <output>...

FLAGS:
    -e, --extract    Extract reads instead of removing them
    -h, --help       Prints help information
    -K, --keep       Keep the working directory and intermediate files
    -V, --version    Prints version information

OPTIONS:
    -i, --input <input>...                                Input filepath(s) (fa, fq, gz, bz)
    -o, --output <output>...                              Output filepath(s) with reads removed or extracted 
    -J, --json <json>                                     Output filepath for summary of depletion/extraction
    -k, --kraken-db <kraken-db>...                        `Kraken2` database directory path(s)
    -t, --kraken-taxa <kraken-taxa>...                    Taxa and sub-taxa (Domain and below) to include
    -d, --kraken-taxa-direct <kraken-taxa-direct>...      Taxa to include directly from reads classified
    -j, --kraken-threads <kraken-threads>                 Threads to use for `Kraken2` [default: 4]
    -m, --minimap2-index <minimap2-index>...              Reference sequence or index file(s) for `minimap2`
    -x, --minimap2-preset <sr|map-ont|map-hifi|map-pb>    `Minimap2` preset configuration [default: sr]
    -n, --minimap2-threads <minimap2-threads>             Threads to use for `minimap2` [default: 4]
    -s, --strobealign-index <strobealign-index>...        Reference sequence (.fa|.fasta) or index (.sti) file(s) for `strobealign`
    -y, --strobealign-preset <map|align>                  `Strobealign` mode configuration [default: align]
    -p, --strobealign-threads <strobealign-threads>       Threads to use for `strobealign` [default: 4]
    -c, --min-cov <min-cov>                               Minimum query alignment coverage to deplete a read [default: 0]
    -l, --min-len <min-len>                               Minimum query alignment length to deplete a read [default: 0]
    -q, --min-mapq <min-mapq>                             Minimum mapping quality to deplete a read [default: 0]
    -O, --output-format <u|b|g|l>                         u: uncompressed; b: Bzip2; g: Gzip; l: Lzma
    -L, --compression-level <1-9>                         Compression level to use if compressing output [default: 6]
    -W, --workdir <workdir>                               Working directory containing intermediary files
```

### Kraken scrubber


```shell
scrubby-scrub-kraken 0.3.0
Deplete or extract reads using outputs from Kraken2

USAGE:
    scrubby scrub-kraken [FLAGS] [OPTIONS] --input <input>... --kraken-reads <kraken-reads> --kraken-report <kraken-report> --output <output>...

FLAGS:
    -e, --extract    Extract reads instead of removing them
    -h, --help       Prints help information
    -V, --version    Prints version information

OPTIONS:
    -L, --compression-level <1-9>                       Compression level to use [default: 6]
    -i, --input <input>...                              Input filepath(s) (fa, fq, gz, bz)
    -J, --json <json>                                   Output filepath for summary of depletion/extraction
    -n, --kraken-name <kraken-name>                     Database name for JSON summary, default is --kraken-reads filestem
    -k, --kraken-reads <kraken-reads>                   Kraken2 classified reads output
    -r, --kraken-report <kraken-report>                 Kraken2 taxonomic report output
    -t, --kraken-taxa <kraken-taxa>...                  Taxa and sub-taxa (Domain and below) to include
    -d, --kraken-taxa-direct <kraken-taxa-direct>...    Taxa to include directly from reads classified
    -o, --output <output>...                            Output filepath(s) with reads removed or extracted
    -O, --output-format <u|b|g|l>                       u: uncompressed; b: Bzip2; g: Gzip; l: Lzma
    -W, --workdir <workdir>                             Working directory for intermediary files
```


### Alignment scrubber


```shell
scrubby-scrub-alignment 0.2.1
Deplete or extract reads using alignments (PAF|SAM|BAM|CRAM)

USAGE:
    scrubby scrub-alignment [FLAGS] [OPTIONS] --alignment <alignment> --input <input>... --output <output>...

FLAGS:
    -e, --extract    Extract reads instead of removing them
    -h, --help       Prints help information
    -V, --version    Prints version information

OPTIONS:
    -a, --alignment <alignment>                    Alignment file (SAM/BAM/CRAM/PAF) or list of read identifiers (TXT)
    -A, --alignment-format <bam|paf|txt|kraken>    bam: SAM/BAM/CRAM alignment; paf: PAF alignment, txt: read identifiers
    -n, --alignment-name <alignment-name>          Alignment name for JSON summary, by default uses --alignment filestem
    -L, --compression-level <1-9>                  Compression level to use [default: 6]
    -i, --input <input>...                         Input filepath(s) (fa, fq, gz, bz)
    -J, --json <json>                              Output filepath for summary of depletion/extraction
    -c, --min-cov <min-cov>                        Minimum query alignment coverage filter [default: 0]
    -l, --min-len <min-len>                        Minimum query alignment length filter [default: 0]
    -q, --min-mapq <min-mapq>                      Minimum mapping quality filter [default: 0]
    -o, --output <output>...                       Output filepath(s) with reads removed or extracted
    -O, --output-format <u|b|g|l>                  u: uncompressed; b: Bzip2; g: Gzip; l: Lzma
    -W, --workdir <workdir
```

## Considerations

### Taxonomic database errors

It should be ensured that the `Kraken2` database correctly specifies taxonomic ranks so that, for example, no further major domain ranks (D) are contained within the domain Eukaryota.

This may be the case in some databases like the [SILVA rRNA](https://benlangmead.github.io/aws-indexes/k2) index which incorrectly specifies Holozoa and Nucletmycea (sub-ranks of domain Eukaryota) as domain (D). Fortunately, it does not appear to be the case for the major [RefSeq databases like PlusPF](https://benlangmead.github.io/aws-indexes/k2).

You can check your database ranks before using `Scrubby`:

```
kraken2-inspect --db silva/ | grep -P "\tD\t"
```

|       |          |         |   |       |             |
|-------|----------|---------|---|-------|-------------|
| 72.87 | 18794250 | 1187248 | D | 3     | Bacteria    |
| 22.42 | 5782938  | 140693  | D | 4     | Eukaryota   |
| 9.28  | 2393894  | 0       | D | 46959 | Holozoa     |
| 2.25  | 580762   | 1017    | D | 47567 | Nucletmycea |
| 4.65  | 1198684  | 44413   | D | 2     | Archaea     |


In this case, Holozoa and Nucletmycea should be added to the `--kraken-taxa` argument, so that all Eukaryota are parsed correctly:

```
scrubby scrub-reads \
  --input R1.fq.gz R2.fq.gz \
  --output S1.fq.gz S2.fq.gz \
  --kraken-db silva/ \
  --kraken-taxa Eukaryota Holozoa Nucletmycea \
  --minimap2-index chm13v2.fasta \
  --min-len 50
```

## Roadmap

* `v0.4.0` - alignment step with `Bowtie2`, `BioConda` deployment for `Linux/OSX`
* `v0.5.0` - metagenome assembly re-alignment depletion with `metaSPAdes` (short reads)
* `v0.6.0` - unit / integration tests

## Dependencies

`Scrubby` wraps or implements the following libraries and tools. If you are using `Scrubby` for publication, please cite:

* [`niffler`](https://github.com/luizirber/niffler)
* [`needletail`](https://github.com/onecodex/needletail)
* [`minimap2`](https://github.com/lh3/minimap2)
* [`strobealign`](https://github.com/ksahlin/strobealign)
* [`Kraken2`](https://github.com/DerrickWood/kraken2)
