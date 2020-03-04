# Transcriptomic analysis of fenitrothion resistance in *Aedes aegypti* from Angola 

Code and data to reproduce the analyses of differential gene expression in pirimiphos-methyl resistant *A. gambiae* and *A. coluzzii* from Côte d'Ivoire. Each of the resistant strains was independently compared to lab colonies from the same species (see below).

## Setup

To re-do the main analyses of the paper, just clone this repository and run the following *R* script:

```bash
git clone git@github.com:xgrau/fenitrothion-aegypti-Angola # or download it
```

## Differential expression

To run the differential expression analysis:

```bash
Rscript XXXXX.R
```

The scripts will:

* load input data,
* perform the differential expression analysis between various groups of samples
* perform functional enrichment analyses of differentially expressed genes
* create figures and tables and various outputs, for each species (`results_de/`).

You can find gene-wise differential expression statistics for each comparison (fold changes, p-values, etc.) in the `results_de/` folder, specifically in the `session.deseq_difexp.ALL.RUN-SRO.csv` tables. There is one table per each comparison:

* `RUN-SRO`: resistant unexposed to susceptible Rockefeller colony. Positive FC = overexpression in RUN.
* `RUN-SNO`: resistant unexposed to susceptible New Orléans colony. Positive FC = overexpression in RUN.
* `RUN-RFE`: resistant unexposed to resistant exposed. Positive FC = overexpression in RUN.

Naming converions for the other output files in `results_de/`:

* files prefixed with `de_all`: results of some analyses of differential expression common to all samples, such as PCAs, sample clustering, heatmaps, etc.
* files prefixed with `de_co`: results from differential expression analyses specific to one comparison. For example, files with `de_co.ALL.RUN-SRO` contain info on the comparison of `RUN` (resistant unexposed) to `SRO` (Susceptible Rockefeller colony) samples.
* files prefixed with `session`: concatenated outputs from the various individual comparisons (and the *R* session).

All necessary input data and gene annotations are included in this package (see **Contents** section).

*R* libraries required to run these analyses are listed below (**Requirements** section).

## Contents

Folders:

* `data_expression/`: sample-specific gene expression data, as produced by *salmon*. It can be loaded into *DESeq2*.
* `data_genome/`: genome annotations (GO, Pfam, gene names, etc.) required for functional enrichments and gene annotation.
* `data_metadata/`: list of samples (12 in total) and sample classification. Sample classification codes are used in *DESeq2* to define specific comparisons between sample groups.
* `helper_scripts`: custom scripts for basic plots (heatmaps, volcano plots) and functional gene enrichment. Sourced by the main scripts.
* `results_de/`: output from differenial expression analysis

## Methods

### Sample list

List of samples:

| Sample code | Category | Resistant? | Exposed? | Colony | Resistant+Exposed? | 
| ------- | --- | --- | ---- | ------ | --- |
| s01_RUN | RUN | Res | Unex | Angola | RUN |
| s02_RUN | RUN | Res | Unex | Angola | RUN |
| s03_RUN | RUN | Res | Unex | Angola | RUN |
| s12_RFE | RFE | Res | Feex | Angola | RFE |
| s13_RFE | RFE | Res | Feex | Angola | RFE |
| s14_RFE | RFE | Res | Feex | Angola | RFE |
| s05_SNO | SNO | Sus | Unex | NewOrl | SUS |
| s06_SNO | SNO | Sus | Unex | NewOrl | SUS |
| s07_SNO | SNO | Sus | Unex | NewOrl | SUS |
| s08_SRO | SRO | Sus | Unex | Rockef | SUS |
| s09_SRO | SRO | Sus | Unex | Rockef | SUS |
| s10_SRO | SRO | Sus | Unex | Rockef | SUS |

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

Reads have been mapped to *Aedes aegypti* transcripts, annotation version L5.1. Downloaded from  [Vectorbase](https://www.vectorbase.org/downloads). Only longest transcript per gene.

### Requirements

Analyses were run on *R* 3.6.1.

The following libraries are required to run the main script:

```R
library(ape)
library(gplots)
library(stringr)
library(plyr)
library(tidyr)
library(topGO)
library(tximport)
library(readr)
library(DESeq2)
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
