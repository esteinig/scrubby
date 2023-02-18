# scrubby

A (t)rusty read scrubber to deplete/extract background taxa using k-mer classifications or alignments.

## Overview

`Scrubby` can deplete/extract reads classified at taxonomic sub-ranks and perform sequential depletion/extraction from multiple databases or alignments. As an example, you can specify a primary (fast) k-mer depletion of all reads classified as _Eukaryota_ (including sub-ranks like _Holozoa_) with `Kraken2`, then follow up with `minimap2` alignment against [`CHM13v2`](https://github.com/marbl/CHM13) to further deplete those pesky human reads.

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

Single or paired-end reads are supported. Compression formats are recognized from extensions of `--input/--output` (`gz|bz|bz2|lzma`). 
Taxa specified for `Kraken2` depletion can be taxonomic identifiers or taxonomic names present in the report file (case sensitive).

```
scrubby scrub-reads \
  --input R1.fq.gz R2.fq.gz \
  --output S1.fq.gz S2.fq.gz \
  --kraken-db minikraken/ \
  --kraken-taxa Eukaryota \
  --minimap2-index chm13v2.fasta \
  --min-len 50
```

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

It should be ensured that the `Kraken2` database correctly specifies taxonomic ranks so that, for example, no further major domain ranks (D) are contained within the domain _Eukaryota_. Minor rank designations are considered to be sub-ranks and depleted within _Eukaryota_ (D1, D2, ...)

This may be the case in some databases like the [SILVA rRNA](https://benlangmead.github.io/aws-indexes/k2) index which incorrectly specifies _Holozoa_ and _Nucletmycea_ (sub-ranks of domain _Eukaryota_) as domain (D). Fortunately, it does not appear to be the case for the major [RefSeq indices](https://benlangmead.github.io/aws-indexes/k2).

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


In this case, _Holozoa_ and _Nucletmycea_ should be added to the `--kraken-taxa` argument, so that all _Eukaryota_ are parsed correctly:

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

* `0.3.0` - secondary alignment depletion with `strobealign` or `Bowtie2`
* `0.4.0` - read extraction and independent sub-commands for alignments and `Kraken2` (from outputs rather than execution)
* `0.5.0` - machine-readable output summaries for workflow integrations, actions for release builds and unit testing

## Dependencies

`Scrubby` wraps or implements the following libraries and tools. If using `Scrubby` in publications, please cite:

* [`niffler`](https://github.com/luizirber/niffler)
* [`needletail`](https://github.com/onecodex/needletail)
* [`minimap2`](https://github.com/lh3/minimap2)
* [`strobealign`](https://github.com/ksahlin/strobealign)
* [`Kraken2`](https://github.com/DerrickWood/kraken2)
