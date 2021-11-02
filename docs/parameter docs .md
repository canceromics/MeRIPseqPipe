# Parameters

## >Input/output options

| Params          | Default | Description |
| --------------- | ------- | ----------- |
| `--designfile`  |         |             |
| `--comparefile` |         |             |
| `--outdir`      |         |             |

## >Data feature options

| Params         | Default | Description |
| -------------- | ------- | ----------- |
| `--single_end` | false   |             |
| `--gzip`       | true    |             |
| `--stranded`   | no      |             |

## >Reference genome options

| Params              | Default | Description |
| ------------------- | ------- | ----------- |
| `--genome`          |         |             |
| `--fasta`           |         |             |
| `--gtf`             |         |             |
| `--star_index`      |         |             |
| `--bwa_index`       |         |             |
| `--hisat2_index`    |         |             |
| `--tophat2_index`   |         |             |
| `--rRNA_fasta`      |         |             |
| `--skip_filterrRNA` | false   |             |

## >Quality control options

| Params                  | Default | Description |
| ----------------------- | ------- | ----------- |
| `--skip_fastp`          | false   |             |
| `--skip_fastqc`         | false   |             |
| `--skip_rseqc`          | false   |             |
| `--skip_qc`             | false   |             |
| `--skip_createbedgraph` | true    |             |

## >Alignment options

| Params          | Default | Description                                        |
| --------------- | ------- | -------------------------------------------------- |
| `--aligners`    | star    | "star" OR "bwa" OR "tophat2" OR "hisat2" OR "none" |
| `--skip_sort`   | false   |                                                    |
| `--mapq_cutoff` |         |                                                    |

## >Expression analysis options

| Params                         | Default | Description |
| ------------------------------ | ------- | ----------- |
| `--featurecounts_feature_type` |         |             |
| `--expression_analysis_mode`   |         |             |

## >Peak calling options

| Params               | Default      | Description                 |
| -------------------- | ------------ | --------------------------- |
| `--peakCalling_mode` | independence | "group" OR "independence"   |
| `--peak_threshold`   | medium       | "low" OR "medium" OR "high" |
| `--skip_peakCalling` | false        |                             |
| `--skip_metpeak`     | false        |                             |
| `--skip_macs2`       | false        |                             |
| `--skip_matk`        | false        |                             |
| `--skip_meyer`       | false        |                             |

## >Peak merging options

| Params              | Default | Description                                        |
| ------------------- | ------- | -------------------------------------------------- |
| `--peakMerged_mode` | rank    | "rank" OR "macs2" OR "MATK" OR "metpeak" OR "mspc" |

## >Methylation analysis options

| Params                        | Default | Description                                                  |
| ----------------------------- | ------- | ------------------------------------------------------------ |
| `--methylation_analysis_mode` | QNB     | "MATK" OR "QNB" OR "Wilcox-test" OR "MeTDiff" OR "edgeR" OR "DESeq2" |
| `--motif_length`              | 5,6,7,8 |                                                              |
| `--skip_annotation`           | false   |                                                              |
| `--skip_m6Aprediction`        | false   |                                                              |
| `--skip_diffpeakCalling`      | false   |                                                              |
| `--skip_motif`                | false   |                                                              |

## >Other options

| Params              | Default                               | Description                                                  |
| ------------------- | ------------------------------------- | ------------------------------------------------------------ |
| `--name`            | false                                 |                                                              |
| `--igenomes_base`   | 's3://ngi-igenomes/igenomes/'         | Directory / URL base for iGenomes references.                |
| `--igenomes_ignore` | false                                 | Do not load the iGenomes reference config.                   |
| `--max_cpus`        | 16                                    | Maximum number of CPUs that can be requested for any single job. |
| `--max_memory`      | 128.GB                                | Maximum amount of memory that can be requested for any single job. |
| `--max_time`        | 240.h                                 | Maximum amount of time that can be requested for any single job. |
| `--help`            | false                                 | Display help text.                                           |
| `--monochrome_logs` | false                                 | Do not use coloured log outputs.                             |
| `--multiqc_config`  | '$baseDir/assets/multiqc_config.yaml' | Custom config file to supply to MultiQC.                     |
| `--tracedir`        | '${params.outdir}/pipeline_info'      | Directory to keep pipeline Nextflow logs and reports.        |

