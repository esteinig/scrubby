# scrubby

A (t)rusty read scrubber to deplete/extract background taxa using k-mer classifications or alignments.

## Overview

**`v0.2.1`**

- [Purpose](#purpose)
- [Install](#install)
- [Usage](#usage)
  - [Input and output](#input-and-output)
    - [General operations](#general-operations)
    - [Read scrubbing](#read-scrubbing)
    - [Kraken2 scrubbing](#kraken2-scrubbing)
    - [Alignment scrubbing](#alignment-scrubbing)
    - [Logging outputs](#logging) 
- [Considerations](#considerations)
  - [Taxonomic sub-rank depletion](#taxonomic-sub-rank-depletion)
  - [Taxonomic database misassignments](#taxonomic-database-misassignments)
- [Command-line arguments](#command-line-arguments)
- [Roadmap](#roadmap)
- [Dependencies](#dependencies)

## Purpose

`Scrubby` can deplete/extract reads classified at taxonomic sub-ranks and perform sequential depletion/extraction from multiple databases or alignments. 

As an example, you can specify a primary (fast) k-mer depletion of all reads classified as Eukaryota (including sub-ranks like Holozoa) with `Kraken2`, then follow up with `minimap2` alignment against [`CHM13v2`](https://github.com/marbl/CHM13) to further deplete those pesky human reads.

`Scrubby` is the mirror counterpart of the excellent [`ReadItAndKeep`](https://github.com/GlobalPathogenAnalysisService/read-it-and-keep) and meant to facilitate safe and thorough background depletion of host and other taxa for metagenomic applications when there is no specific target genome. However, you can also use `Scrubby` for target retention by specifying the target genome/taxon and `--extract` flag to retain reads aligned/classified by one or multiple methods (v0.3.0) .

This is a preliminary release, use at your own peril :skull:

## Install

BioConda:

```
conda install -c conda-forge -c bioconda scrubby
```

Cargo:

```
cargo install scrubby
```

## Usage

Expanded command-line argument descriptions are available with the `--help` flag.


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

### Input and Output

#### General operations

Single or paired-end reads are supported. Compression formats are recognized from extensions of `--input/--output` (`gz|bz|bz2|xz`). 
Taxa specified for `Kraken2` depletion can be taxonomic identifiers or taxonomic names present in the report file (case sensitive, see below).
Alignment filters can can be specified (using query alignment length, coverage or mapping quality as in `ReadItAndKeep`). 

A JSON formatted read depletion/extraction summary can be written to file (`--json file.json`) or stdout (`--json -`). The schema contains a summary 
array for each database or reference provided in the order in which reads were depleted/extracted in the scrubbing pipeline. When individually depleting/extracting
`Kraken2` (`scrubby scrub-kraken`) or alignments (`scrubby scrub-alignments`), this is always an index of `0` and the name of the `Kraken2` reads file
or the name of the alignment.

```json
{
  "version": "0.2.1",
  "schema_version": "0.1.0",
  "summary": [
    {
      "index": 0,
      "name": "rrna",
      "total": 1000000,
      "depleted": 66784,
      "retained": 933216,
      "extracted": 0,
      "files": [
        {
          "total": 500000,
          "depleted": 33392,
          "retained": 466608,
          "extracted": 0,
          "input_file": "/path/to/test_r1.fq",
          "output_file": "/path/to/tmp/workdir/0-rrna_1.fq"
        },
        {
          "total": 500000,
          "depleted": 33392,
          "retained": 466608,
          "extracted": 0,
          "input_file": "/path/to/test_r2.fq",
          "output_file": "/path/to/tmp/workdir/0-rrna_2.fq"
        }
      ]
    }
  ]
}
```

#### Read scrubbing

`Scrubby` primarily depletes/extracts using sequential operation of k-mer and alignment methods. This will call `Kraken2` and aligners under the hood,
creating intermediary files in the `-W/--wordir` which can be retained (`-K/--keep`). By default the working directory is created with
a time-stamp (`Scrubby_{YYYYMMDDTHHMMSS}`).


```
scrubby scrub-reads \
  --input R1.fq.gz R2.fq.gz \
  --output S1.fq.gz S2.fq.gz \
  --kraken-db minikraken/ \
  --kraken-taxa Metazoa \
  --minimap2-index chm13v2.fasta \
  --min-len 50
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

#### Logging

Logs are output to stderr:

![output](https://user-images.githubusercontent.com/12873366/219830790-03deeb50-40de-4587-bff2-0111cc620300.png)


## Considerations

### Taxonomic sub-rank depletion 

When using `Kraken2`, all reads classified within a particular taxonomic rank **except those above Domain** can be depleted/extracted. If you only want to deplete/extract a specific taxon or want to deplete/extract a rank above Domain use the `--kraken-taxa-direct` argument (for example "Unclassified" and "Cellular Organisms") :

```
scrubby scrub-reads \
  --input R1.fq.gz R2.fq.gz \
  --output S1.fq.gz S2.fq.gz \
  --kraken-db minikraken/ \
  --kraken-taxa Eukaryota \
  --kraken-taxa-direct Unclassified 131567
```

### Taxonomic database misassignments

It should be ensured that the `Kraken2` database correctly specifies taxonomic ranks so that, for example, no further major domain ranks (D) are contained within the domain Eukaryota. Minor rank designations are considered to be sub-ranks and depleted within Eukaryota (D1, D2, ...)

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

## Command-line arguments

```shell
scrubby-scrub-reads 0.2.1
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
    -k, --kraken-db <kraken-db>...                        Kraken2 database directory path(s)
    -t, --kraken-taxa <kraken-taxa>...                    Taxa and sub-taxa (Domain and below) to include
    -d, --kraken-taxa-direct <kraken-taxa-direct>...      Taxa to include directly from reads classified
    -j, --kraken-threads <kraken-threads>                 Threads to use for Kraken2 [default: 4]
    -m, --minimap2-index <minimap2-index>...              Reference sequence or index file(s) for `minimap2`
    -x, --minimap2-preset <sr|map-ont|map-hifi|map-pb>    Minimap2 preset configuration [default: sr]
    -n, --minimap2-threads <kraken-threads>               Threads to use for minimap2 [default: 4]
    -c, --min-cov <min-cov>                               Minimum query alignment coverage to deplete a read [default: 0]
    -l, --min-len <min-len>                               Minimum query alignment length to deplete a read [default: 0]
    -q, --min-mapq <min-mapq>                             Minimum mapping quality to deplete a read [default: 0]
    -O, --output-format <u|b|g|l>                         u: uncompressed; b: Bzip2; g: Gzip; l: Lzma
    -L, --compression-level <1-9>                         Compression level to use if compressing output [default: 6]
    -W, --workdir <workdir>                               Working directory containing intermediary files
```

## Roadmap

* `0.3.0` - secondary alignment depletion with `strobealign` or `Bowtie2`, another attempt to compile `rust-htslib` fail for Bioconda OSX
* `0.4.0` - read extraction and independent sub-commands for alignments and `Kraken2` (from outputs rather than execution)
* `0.5.0` - machine-readable output summaries for workflow integrations, actions for release builds and unit testing

## Dependencies

`Scrubby` wraps or implements the following libraries and tools. If you are using `Scrubby` for publication, please cite:

* [`niffler`](https://github.com/luizirber/niffler)
* [`needletail`](https://github.com/onecodex/needletail)
* [`minimap2`](https://github.com/lh3/minimap2)
* [`strobealign`](https://github.com/ksahlin/strobealign)
* [`Kraken2`](https://github.com/DerrickWood/kraken2)
