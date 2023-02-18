# scrubby

A (t)rusty read scrubber to deplete/extract background reads using taxonomic classifications or filtered alignments.

## Overview

`Scrubby` can be used to deplete/extract reads classified at sub-taxonomic levels and perform sequential depletion/extraction from multiple databases or alignments. Initial release supports `minimap2` and `Kraken2`. As an example, you can specify a primary (fast) k-mer depletion of all reads classified as `Eukaryota` (including sub-taxonomic levels) with `Kraken2`, then follow up with `minimap2` alignment against [`CHM13v2`](https://github.com/marbl/CHM13) to further deplete those pesky human reads.

This is a preliminary release, use at your own peril :skull:

## Install

BioConda:

```
conda install -c conda-forge -c bioconda scrubby
```

Cargo package manager:

```
cargo install scrubby
```

## Usage

Expanded command-line argument descriptions:

```
scrubby scrub-reads --help
```

Single or paired-end reads are supported. Compression formats are recognized from file extensions of `--input/--output`

```
scrubby scrub-reads \
  --input R1.fq.gz R2.fq.gz \
  --output S1.fq.gz S2.fq.gz \
  --kraken-db minikraken/ \
  --kraken-taxa Eukaryota \
  --minimap2-index chm13v2.fasta \
  --min-len 50
```

### Taxonomic sub-level depletion

It should be ensured that the `Kraken2` database correctly specifies taxonomic levels so that, for example, no further major domain levels (D) are contained within Eukaryota. Minor level designations are considered to be sub-levels and depleted within Eukaryota (D1, D2, ...)

This may be the case in some databases like the [SILVA rRNA](https://benlangmead.github.io/aws-indexes/k2) index which incorrectly specifies `Holozoa` and `Nucletmycea` (sub-levels of domain Eukaryota) as domain (D). 

You can check your database levels before using `Scrubby`:

```
kraken2-inspect --db silva/ | grep -P "\tD\t"
```

In this case, `Holozoa` and `Nucletmycea` should be added to the `--kraken-taxa` argument, so that all `Eukaryota` are parsed correctly:

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
* `0.5.0` - machine-readable output summaries for workflow integrations

## Dependencies

`Scrubby` wraps or implements the following libraries and tools - please cite if used in publications:

* [`niffler`](https://github.com/luizirber/niffler)
* [`needletail`](https://github.com/onecodex/needletail)
* [`minimap2`](https://github.com/lh3/minimap2)
* [`strobealign`](https://github.com/ksahlin/strobealign)
* [`Kraken2`](https://github.com/DerrickWood/kraken2)
