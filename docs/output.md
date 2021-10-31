# MeRIPseqPipe: Output

This document describes the output produced by the pipeline. Most of the plots are taken from the MultiQC report, which summarises results at the end of the pipeline.

<!-- TODO nf-core: Write this documentation describing your workflow's output -->

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/)
and processes data using the following steps:

- [MeRIPseqPipe: Output](#MeRIPseqPipe-output)
  - [Pipeline overview](#pipeline-overview)
  - [Quality Control](#quality-control)
    - [Fastp](#fastp)
    - [FastQC](#fastqc)
    - [RSeQC](#rseqc)
    - [MultiQC](#multiqc)
  - [Align reasults](#align-reasults)
    - [STAR](#star)
    - [BWA](#bwa)
    - [TopHat2](#tophat2)
    - [HISAT2](#hisat2)
  - [SAMtools](#samtools)
  - [PeakCalling](#peakcalling)
    - [MeTPeak](#metpeak)
    - [MATK](#matk)
    - [Meyer](#meyer)
    - [MACS2](#macs2)
  - [PeakMerged](#peakmerged)
    - [RobustRankAggreg](#robustrankaggreg)
    - [MSPC](#mspc)
    - [BEDtools](#bedtools)
  - [M6A sites prediction](#m6a-sites-prediction)
    - [MATK](#matk)
  - [Differtial Methylation Analysis](#differtial-methylation-analysis)
    - [QNB](#qnb)
    - [MATK](#matk)
    - [DESeq2_DM](#deseq2dm)
    - [edgeR_DM](#edgerdm)
  - [Differtial Expression Analysis](#differtial-expression-analysis)
    - [featureCounts](#featureCounts)
    - [DESeq2_DE](#deseq2de)
    - [edgeR_DE](#edgerde)
  - [Reports](#reports)

Several R packages for downstream analysis.

## Quality Control

**Output directory: `results/QC/`**

### Fastp

[Fastp](https://github.com/OpenGene/fastp)

**Output directory: `results/QC/fastp`**

### FastQC

[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) gives general quality metrics about your reads. It provides information about the quality score distribution across your reads, the per base sequence content (%T/A/G/C). You get information about adapter contamination and other overrepresented sequences.

For further reading and documentation see the [FastQC help](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/).

> **NB:** The FastQC plots displayed in the MultiQC report shows _untrimmed_ reads. They may contain adapter sequence and potentially regions with low quality. To see how your reads look after trimming, look at the FastQC reports in the `trim_galore` directory.

**Output directory: `results/QC/fastqc`**

- `sample_fastqc.html`
  - FastQC report, containing quality metrics for your untrimmed raw fastq files
- `zips/sample_fastqc.zip`
  - zip file containing the FastQC report, tab-delimited data file and plot images

### RSeQC

[RSeQC](http://rseqc.sourceforge.net/)

**Output directory: `results/QC/RSeQC`**

### MultiQC

[MultiQC](http://multiqc.info) is a visualisation tool that generates a single HTML report summarising all samples in your project. Most of the pipeline QC results are visualised in the report and further statistics are available in within the report data directory.

The pipeline has special steps which allow the software versions used to be reported in the MultiQC output for future traceability.

**Output directory: `results/QC/multiqc`**

- `Project_multiqc_report.html`
  - MultiQC report - a standalone HTML file that can be viewed in your web browser
- `Project_multiqc_data/`
  - Directory containing parsed statistics from the different tools used in the pipeline

For more information about how to use MultiQC reports, see [http://multiqc.info](http://multiqc.info)

## Align reasults

**Output directory: `results/QC/multiqc`**

### STAR

[STAR](https://github.com/alexdobin/STAR)

### BWA

[BWA](https://github.com/lh3/bwa)

### TopHat2

[TopHat2](https://ccb.jhu.edu/software/tophat/)

### HISAT2

[HISAT2](https://ccb.jhu.edu/software/hisat2/)

## SAMtools

[SAMtools](http://www.htslib.org/)

## PeakCalling

### MeTPeak

[MeTPeak](https://github.com/compgenomics/MeTPeak)

### MATK

[MATK](http://matk.renlab.org)

### Meyer

[meyer](http://matk.renlab.org)

### MACS2

[MACS2](https://github.com/taoliu/MACS)

## PeakMerged

### RobustRankAggreg

[RobustRankAggreg](https://cran.r-project.org/web/packages/RobustRankAggreg/index.html)

### MSPC

[MSPC]

### BEDtools

[BEDtools](https://bedtools.readthedocs.io/en/latest/index.html)

## M6A sites prediction

[MATK](http://matk.renlab.org)

## Differtial Methylation Analysis

### QNB

[QNB](https://cran.r-project.org/src/contrib/Archive/QNB/)

### MATK
[MATK](http://matk.renlab.org)

### DESeq2_DM

[DESeq2](http://bioconductor.org/packages/DESeq2/)

### edgeR_DM

[edgeR](http://bioconductor.org/packages/edgeR/)

## Differtial Expression Analysis

### featureCounts

[featureCounts](http://subread.sourceforge.net)

### DESeq2_DE

[DESeq2](http://bioconductor.org/packages/DESeq2/)

### edgeR_DE

[edgeR](http://bioconductor.org/packages/edgeR/)

## Reports
