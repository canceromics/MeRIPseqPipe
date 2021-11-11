# Parameters

- [Default parameters](#default-parameters)
- [Input/output options](#inputoutput-options)
- [Data feature options](#data-feature-options)
- [Reference genome options](#feference-genome-options)
- [Quality control options](#quality-control-options)
- [Alignment options](#alignment-options)
- [Expression analysis options](#expression-analysis-options)
- [Peak calling options](#peak-calling-options)
- [Peak merging options](#peak-merging-options)
- [Methylation analysis options](#methylation-analysis-options)
- [Plot options](#plot-options)
- [Other options](#other-options)

## >Default parameters

| Tools or packages or methods                   | Version | Parameters (taking Single-end data as example)               |
| ---------------------------------------------- | ------- | ------------------------------------------------------------ |
| **Raw data quality control and preprocessing** |         |                                                              |
| fastp                                          | 0.19.7  | default                                                      |
| FastQC                                         | 0.11.8  | default                                                      |
| **rRNA filtering**                             |         |                                                              |
| HISAT2                                         | 2.1.0   | `hisat2 --summary-file ${sample_name}_rRNA_summary.txt --no-spliced-alignment --no-softclip --norc --no-unal -p ${task.cpus} --dta --un-gz ${sample_id}.fastq.gz -x genome_index -U read1.fastq.gz` |
| **Reads alignment**                            |         |                                                              |
| STAR                                           | 2.6.1b  | `STAR --runThreadN ${task.cpus} --readFilesCommand zcat --twopassMode Basic --genomeDir genome_index --readFilesIn read1.fastq.gz --outSAMtype BAM Unsorted --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterIntronMotifs RemoveNoncanonical --outFilterMultimapNmax 20 --alignIntronMin 20 --alignIntronMax 100000 --alignMatesGapMax 1000000 --outFileNamePrefix ${sample_name} > ${sample_name}_log.txt` |
| HISAT2                                         | 2.1.0   | `hisat2 --summary-file ${sample_name}_hisat2_summary.txt -p ${task.cpus} --dta -x genome_index -U read1.fastq.gz` |
| BWA                                            | 0.7.17  | `bwa aln -t ${task.cpus} -f ${read1.baseName}.sai genome_index read1.fastq.gz <br />bwa samse genome_index ${read1.baseName}.sai read1.fastq.gz` |
| **Quality control after alignment**            |         |                                                              |
| SAMtools                                       | 1.9     | default for converting sam to bam format followed by sorting and indexing, 20 for mapping quality filtering |
| RseQC                                          | 2.6.4   | infer_experiment.py<br />junction_annotation.py<br />bam_stat.py<br />junction_saturation.py<br />inner_distance.py<br />read_distribution.py<br />read_duplication.py |
| **Expression quantification**                  |         |                                                              |
| featureCounts                                  | 2.0.0   | default, `-Q 0 -s 0`                                         |
| **Differential expression analysis**           |         |                                                              |
| DESeq2                                         | 1.22.1  | default                                                      |
| edgeR                                          | 3.26.0  | default                                                      |
| **Peak calling**                               |         |                                                              |
| MeTPeak                                        | 1.0.0   | default                                                      |
| MACS2                                          | 2.1.2   | default, `-q 0.01 --keep-dup 5 -f BAM --nomodel`             |
| MATK                                           | 1.0     | default, `-q 0.05`                                           |
| Meyer                                          |         | python scripts                                               |
| **Peak merging**                               |         |                                                              |
| RobustRankAggreg                               | 1.1     | R scripts                                                    |
| MSPC                                           | 5.4.0   | default, `-r bio -s 1E-4 -w 1E-2 / -r tec -s 1E-2 -w 1E-1`   |
| BEDTools                                       | 2.27.1  | intersectBed                                                 |
| Peak annotated                                 |         | Perl scripts                                                 |
| **motif searching**                            |         |                                                              |
| HOMER                                          | 4.9.1   | defalut, `-len 5,6,7,8 -S 10 -rna -dumpFasta`                |
| **methylation quantification**                 |         |                                                              |
| QNB                                            | 1.1.11  | mode="blind"                                                 |
| MATK                                           | 1.0     | -quantification                                              |
| RPKM method                                    | 2.27.1  | (ip_rpkm+1)/(input_rpkm+1)                                   |
| **Differential methylation analysis**          |         |                                                              |
| QNB                                            | 1.1.11  | mode="auto"                                                  |
| MATK                                           | 1.0     | -diff                                                        |
| DESeq2(GLM)                                    | 1.22.1  | R scripts, diff.log2fc = peak.log2fc - gene.log2fc           |
| edgeR(GLM)                                     | 3.26.0  | R scripts, diff.log2fc = peak.log2fc - gene.log2fc           |
| Wilcox-test                                    |         | R scripts                                                    |

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
| `--single_end` | false   | "false" for paired_end sequencing;<br />"true" for single_end sequencing. |
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

| Params           | Default | Description                                                  |
| ---------------- | ------- | ------------------------------------------------------------ |
| `--aligners`     | star    | Specifies the alignment algorithm to use.<br />Available option:<br />"star" OR "bwa" OR "tophat2" OR "hisat2" OR "none" |
| `--skip_sort`    | false   | Skip the process of sorting BAM files.                       |
| `--mapq_cutoff ` | 20      | Range from 0 to 255, "255" means the pipeline will only keep the unique mapping reads, "0" means keeping all reads. |

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

## >Plot options

| Params           | Default | Description                                                  |
| ---------------- | ------- | ------------------------------------------------------------ |
| --delfc          | 0.58    | The cutoff of log2 fold change in differential gene filtering. |
| --dmlfc          | 0.58    | The cutoff of log2 fold change in differential peak filtering. |
| --cluster_method | single  | The method of Hierarchical clustering.<br />Available options:<br />"ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid" |
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
