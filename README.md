# scrubby

A (t)rusty read scrubber to deplete/extract background taxa using k-mer classifications or alignments.

## Overview

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

Expanded command-line argument descriptions:

```
scrubby scrub-reads --help
```

Single or paired-end reads are supported. Compression formats are recognized from extensions of `--input/--output` (`gz|bz|bz2|xz`). 
Taxa specified for `Kraken2` depletion can be taxonomic identifiers or taxonomic names present in the report file (case sensitive).

```
scrubby scrub-reads \
  --input R1.fq.gz R2.fq.gz \
  --output S1.fq.gz S2.fq.gz \
  --kraken-db minikraken/ \
  --kraken-taxa Metazoa \
  --minimap2-index chm13v2.fasta \
  --min-len 50
```

A JSON formatted read depletion/extraction summary can be written to file (`--json file.json`) or stdout (`--json -`). The schema contains a summary array for each database or reference provided in the order in which reads were depleted/extract, in this example for paired-end data:

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

Logs are output to `stderr`:

![output](https://user-images.githubusercontent.com/12873366/219830790-03deeb50-40de-4587-bff2-0111cc620300.png)


### Optional methods and order of execution

You can skip methods by leaving out `--kraken-db` or `--minimap2-index` arguments, for example:

```
scrubby scrub-reads \
  --input R1.fq.gz R2.fq.gz \
  --output S1.fq.gz S2.fq.gz \
  --kraken-db minikraken/ \
  --kraken-taxa Eukaryota
```

Order of databases/references used for depletion is the same as the order of the command-line arguments. Order of methods of depletion is currently always: `Kraken2` --> `minimap2` --> `bowtie2` --> `strobealign` (unless a method is not specified).

### Taxonomic sub-ranks in `Kraken2`

`Scrubby` enables depletion of all reads classified within a particular taxonomic rank **except those above Domain**. If you only want to deplete a specific taxon or want to deplete a rank above Domain use the `--kraken-taxa-direct` argument (for example "Unclassified" and "Cellular Organisms") :

```
scrubby scrub-reads \
  --input R1.fq.gz R2.fq.gz \
  --output S1.fq.gz S2.fq.gz \
  --kraken-db minikraken/ \
  --kraken-taxa Eukaryota \
  --kraken-taxa-direct Unclassified 131567
```

### Taxonomic errors in `Kraken2`

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

## Roadmap

* `0.3.0` - secondary alignment depletion with `strobealign` or `Bowtie2`, another attempt to compile `rust-htslib` fail for Bioconda OSX
* `0.4.0` - read extraction and independent sub-commands for alignments and `Kraken2` (from outputs rather than execution)
* `0.5.0` - machine-readable output summaries for workflow integrations, actions for release builds and unit testing

## Dependencies

`Scrubby` wraps or implements the following libraries and tools. If using `Scrubby` in publications, please cite:

* [`niffler`](https://github.com/luizirber/niffler)
* [`needletail`](https://github.com/onecodex/needletail)
* [`minimap2`](https://github.com/lh3/minimap2)
* [`strobealign`](https://github.com/ksahlin/strobealign)
* [`Kraken2`](https://github.com/DerrickWood/kraken2)
