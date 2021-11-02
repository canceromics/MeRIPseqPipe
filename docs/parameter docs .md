# Parameters

## >Input/output options

| Params          | Default                                                   | Description                                                  |
| --------------- | --------------------------------------------------------- | ------------------------------------------------------------ |
| `--designfile`  | “$baseDir/test_datasets/inputfiles/designfile_paired.tsv” | A tab-separated table describing the path to the Input files and IP files (FASTQ or BAM format) and the grouping information (like “control” or “treated”). It's recommended to edit by Excel and save as .tsv suffix file. |
| `--comparefile` | "$baseDir/test_datasets/inputfiles/comparefile.txt"       | A text file that specifies the comparison between two different groups.<br />Other options:<br />"false" for projects without differential analysis;<br />"two_groups" for only two groups in the designfile. |
| `--outdir`      | "$baseDir/results"                                        | Path to the output directory where the results will be saved. |
| `-profile`      | docker                                                    | Configuration files can contain the definition of one or more profiles.<br />Available options:<br />conda, docker, test |

## >Data feature options

| Params         | Default | Description                                                  |
| -------------- | ------- | ------------------------------------------------------------ |
| `--single_end` | false   | "false" for paired_end sequencing;
"true" for single_end sequencing. |
| `--gzip`       | true    | Whether FASTQ files are be gzipped.                          |
| `--stranded`   | no      | "yes" means stranded sequencing data;<br />"no" means unstranded sequncing data;<br />"reverse" means reversely stranded sequencing data. |

## >Reference genome options

| Params              | Default                | Description                                       |
| ------------------- | ---------------------- | ------------------------------------------------- |
| `--genome`          | false                  | Name of iGenomes reference, eg. "GRCh38"          |
| `--fasta`           | "path/to/xxx.fa"       | Path to FASTA genome file.                        |
| `--gtf`             | "path/to/xxx.gtf"      | Path to GTF annotation file.                      |
| `--star_index`      | "path/to/STARIndex"    | Path to directory for pre-built STAR index.       |
| `--bwa_index`       | "path/to/BWAIndex"     | Path to directory for pre-built BWA index.        |
| `--hisat2_index`    | "path/to/HISTA2Index"  | Path to directory for pre-built HISAT2 index.     |
| `--tophat2_index`   | "Path/to/Tophat2Index" | Path to directory for pre-built Tophat2 index.    |
| `--rRNA_fasta`      | "Path/to/rRNA.fa"      | Path to directory for rRNA FASTA genome file.     |
| `--skip_filterrRNA` | false                  | Skip the process of rRNA cleaning before mapping. |

## >Quality control options

| Params                  | Default | Description                                                  |
| ----------------------- | ------- | ------------------------------------------------------------ |
| `--skip_fastp`          | false   | Skip fastp and data preprocessing.                           |
| `--skip_fastqc`         | false   | Skip FastQC.                                                 |
| `--skip_rseqc`          | false   | Skip RseQC.                                                  |
| `--skip_qc`             | false   | Skip all QC steps except for MultiQC.                        |
| `--skip_createbedgraph` | true    | Skip the process of generating bedgraph format files and creatIGVjs |

## >Alignment options

| Params          | Default | Description                                                  |
| --------------- | ------- | ------------------------------------------------------------ |
| `--aligners`    | star    | Specifies the alignment algorithm to use.<br />Available option:<br />"star" OR "bwa" OR "tophat2" OR "hisat2" OR "none" |
| `--skip_sort`   | false   | Skip the process of sorting BAM files.                       |
| `--mapq_cutoff` | 20      | Range from 0 to 255, "255" means the pipeline will only keep the unique mapping reads, "0" means keeping all reads. |

## >Expression analysis options

| Params                       | Default | Description                                                  |
| ---------------------------- | ------- | ------------------------------------------------------------ |
| `--featurecount_minMQS`      | 0       | Integer giving the minimum mapping quality score a read must satisfy in order to be counted. For paired-end reads, at least one end should satisfy this criteria. |
| `--expression_analysis_mode` | DESeq2  | Specifies the tool to do the differential expression analysis.<br />Available options:<br />"DESeq2" OR "edgeR" OR "none" |

## >Peak calling options

| Params               | Default      | Description                                                  |
| -------------------- | ------------ | ------------------------------------------------------------ |
| `--peakCalling_mode` | independence | "group" takes all biological repetitions as a sample,only support MeTPeak and MATK; "independence" takes every sample single. |
| `--peak_threshold`   | medium       | The threshold of peak calling.<br />Available options:<br />"low" OR "medium" OR "high" |
| `--skip_peakCalling` | false        | Skip the process of peak calling.                            |
| `--skip_metpeak`     | false        | Skip MeTPeak.                                                |
| `--skip_macs2`       | false        | Skip MACS2.                                                  |
| `--skip_matk`        | false        | Skip MATK.                                                   |
| `--skip_meyer`       | false        | Skip Meyer.                                                  |

## >Peak merging options

| Params              | Default | Description                                                  |
| ------------------- | ------- | ------------------------------------------------------------ |
| `--peakMerged_mode` | rank    | Specifies  the tool to merge all peaks.<br />Available options:<br />"rank" OR "macs2": merging all peaks from different peak callers.<br />"MATK" OR "metpeak" OR "mspc" OR "meyer": taking the results of this tool as a baseline to merge all peaks. |

## >Methylation analysis options

| Params                        | Default | Description                                                  |
| ----------------------------- | ------- | ------------------------------------------------------------ |
| `--methylation_analysis_mode` | QNB     | Specifies the tool to do the differential methylation analysis.<br />Available options:<br />"MATK" OR "QNB" OR "Wilcox-test" OR "edgeR" OR "DESeq2" |
| `--motif_length`              | 5,6,7,8 | Length for HOMER motif searching.                            |
| `--skip_annotation`           | false   | Skip the process of peak annotation.                         |
| `--skip_m6Aprediction`        | false   | Skip the process of m6A site predicition.                    |
| `--skip_diffpeakCalling`      | false   | Skip all differential methylation analysis.                  |
| `--skip_motif`                | false   | Skip the process of motif searching.                         |

## >Other options

| Params              | Default                               | Description                                                  |
| ------------------- | ------------------------------------- | ------------------------------------------------------------ |
| `--name`            | false                                 | Name for the pipeline run. If not sepcified, nextflow will automatically generate a random mnemonic. |
| `--igenomes_base`   | “s3://ngi-igenomes/igenomes/”         | Directory / URL base for iGenomes references.                |
| `--igenomes_ignore` | false                                 | Do not load the iGenomes reference config.                   |
| `--max_cpus`        | 16                                    | Maximum number of CPUs that can be requested for any single job. |
| `--max_memory`      | 128.GB                                | Maximum amount of memory that can be requested for any single job. |
| `--max_time`        | 240.h                                 | Maximum amount of time that can be requested for any single job. |
| `--help`            | false                                 | Display help text.                                           |
| `--monochrome_logs` | false                                 | Do not use coloured log outputs.                             |
| `--multiqc_config`  | “$baseDir/assets/multiqc_config.yaml” | Custom config file to supply to MultiQC.                     |
| `--tracedir`        | “${params.outdir}/pipeline_info”      | Directory to keep pipeline Nextflow logs and reports.        |
