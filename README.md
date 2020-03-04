# Transcriptomic analyisis of pirmiphos-methyl resistance in *A. gambiae* and *A. coluzzii*

Code and data to reproduce the analyses of differential gene expression in pirimiphos-methyl resistant *A. gambiae* and *A. coluzzii* from Côte d'Ivoire. Each of the resistant strains was independently compared to lab colonies from the same species (see below).

## Setup

To re-do the main analyses of the paper, just clone this repository and run the following `R` script:

```bash
git clone git@github.com:xgrau/pm-resistance-gamcol # or download it
```

## Contents

Folders:

* `data_expression/`: sample-specific gene expression data, as produced by `salmon`. It's formatted for easy loading into `DESeq2`
* `data_genome/`: genome annotations (GO, Pfam, gene names, etc.) required for functional enrichments and gene annotation.
* `data_metadata/`: list of samples (12 in total) and sample classification. Sample classification codes are used in `DESeq2` to define specific comparisons between sample groups.
* `helper_scripts`: custom scripts for basic plots (heatmaps, volcano plots) and functional gene enrichment. Sourced by the main scripts.
* `results*/`: various folders with results, see below.

## Analyses

### Differential expression analyses

To run the differential expression analyses for *gambiae* and *coluzzii* samples:

```bash
# coluzzii differential expression analysis:
Rscript s01_dexp_deseq2_col_2020-01-14.R
# gambiae differential expression analysis:
Rscript s01_dexp_deseq2_gam_2020-01-14.R
```

The scripts will:

* load input data,
* perform the differential expression analysis between various groups of samples from a given species (susceptible colonies or `Col`, resistant controls or `RCo`, and resistant & exposed or `REx`; see below)
* Two types of comparisons: incremental expression in resistant (`Col<RCo<REx`) and consensus over-expression in resistant (`RCo>Col && REx>Col`).
* perform functional enrichment analyses of differentially expressed genes
* create figures and tables and various outputs, for each species (`results_de_gam/` and `results_de_col` folders).

All necessary input data and gene annotations are included in this package (see **Contents** section below).

`R` libraries required to run these analyses are listed below (**Requirements** section).

### Differential splicing analyses

TODO

```bash
# coluzzii differential expression analysis:
Rscript s02_differential_isoforms_gam_2020-02-06.R
# gambiae differential expression analysis:
Rscript s02_differential_isoforms_gam_2020-02-06.R
```

### Genetic variation in the transcriptome

See `README.md` in `results_variation/` folder.

## Methods

### Sample list

For *gambiae* samples

| Sample code | Category | Species | Susceptible? | Exposed? |
|----- | -- | -- | -- | -- |
| Ki11 | Col | gam | Sus | Une |
| Ki2 | Col | gam | Sus | Une |
| Ki5 | Col | gam | Sus | Une |
| Ki7 | Col | gam | Sus | Une |
| SC10 | RCo | gam | Res | Une |
| SC6 | RCo | gam | Res | Une |
| SC7 | RCo | gam | Res | Une |
| SC9 | RCo | gam | Res | Une |
| SS11 | REx | gam | Res | Exp |
| SS5 | REx | gam | Res | Exp |
| SS7 | REx | gam | Res | Exp |
| SS9 | REx | gam | Res | Exp |

*coluzzii* samples:

| Sample code | Category | Species | Susceptible? | Exposed? |
|----- | -- | -- | -- | -- |
| Ng1 | Col | col | Sus | Une |
| Ng11 | Col | col | Sus | Une |
| Ng4 | Col | col | Sus | Une |
| Ng7 | Col | col | Sus | Une |
| MC6 | RCo | col | Res | Une |
| MC7 | RCo | col | Res | Une |
| MC8 | RCo | col | Res | Une |
| MC9 | RCo | col | Res | Une |
| MS1 | REx | col | Res | Exp |
| MS3 | REx | col | Res | Exp |
| MS5 | REx | col | Res | Exp |
| MS7 | REx | col | Res | Exp |

### Read mapping

The process of mapping of reads to the predicted has been performed with `Salmon` v0.10.2, as specified in the Methods section of the paper.

If you want to repeat the analysis yourself, these are the relevant commands (adjusting filenames accordingly):

```bash
# salmon index
salmon index -t transcripts.fasta -i transcripts.salmon.index --type quasi -k 31 1> log_index.log 2> &1
# run this command for each sample separately (A, B, etc.)
salmon quant -i Anogam_long.cds_mcherry.salmon.index -l A -p 10 -1 sampleA_1.fastq.gz -2 sampleA_2.fastq.gz -o sampleA_salmon_out 1> log_quant.log 2> &1
# output in sampleA_salmon_out, sampleB_salmon_out, etc. folders can be loaded to DESeq2 using the main script
```

The read files (fastq format) are not provided in this repository, but they can be found in the ENA public repository under the `XXX` accession number.

Reads have been mapped to *Anopheles gambiae* transcripts, annotation AgamP4.9. Downlodaed from  [Vectorbase](https://www.vectorbase.org/downloads).

## Requirements

Analyses for the paper were run on `R` 3.6.1.

The following libraries are required to run the main script:

```R
library(DESeq2)
library(tximport)
library(ape)
library(gplots)
library(stringr)
library(plyr)
library(tidyr)
library(topGO)
library(readr)
library(pheatmap)
library(cluster)
library(VennDiagram)
```

If you reuse this code in any way, don't forget to cite the original papers for each of these libraries, too. It's easy and makes everyone happy. For example:

```R
> citation("DESeq2")

  Love, M.I., Huber, W., Anders, S. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2 Genome Biology 15(12):550
  (2014)

A BibTeX entry for LaTeX users is

  @Article{,
    title = {Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2},
    author = {Michael I. Love and Wolfgang Huber and Simon Anders},
    year = {2014},
    journal = {Genome Biology},
    doi = {10.1186/s13059-014-0550-8},
    volume = {15},
    issue = {12},
    pages = {550},
  }
```
