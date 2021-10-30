#!/usr/bin/env nextflow

/*
========================================================================================
                         MeRIPseqPipe
========================================================================================
MeRIPseqPipe：An integrated analysis pipeline for MeRIPseq data based on Nextflow.
 #### Homepage / Documentation
 https://github.com/canceromics/MeRIPseqPipe
----------------------------------------------------------------------------------------
*/

/*
 * Authors:
 * Qi Zhao <zhaoqi@sysucc.org.cn>;
 * Bao Xiaoqiong <baoxq@sysucc.org.cn>
 * Zhu Kaiyu <zhuky5@mail2.sysu.edu.cn>:
 */

/* 
requirement:
- FastQC/fastp
- STAR/HISAT2/TopHat2/BWA
- SAMtools/RseQC
- MeTPeak/MACS2/MATK/Meyer
- RobustRankAggreg/MSPC/BEDTools
- featureCounts/DESeq2/EdgeR
- QNB
- HMOER
- igvtools
*/

def helpMessage() {
    log.info"""
    Usage:

    The typical command for running the pipeline is as follows:

    nextflow run path/to/MeRIPseqPipe --designfile 'path/to/designfile/designfile.tsv' --comparefile 'path/to/comparefile/compare.txt' -profile <docker/test/conda> -resume
    
    OR, 

    nextflow path/to/MeRIPseqPipe/main.nf -c path/to/MeRIPseqPipe/nextflow.config --designfile 'path/to/designfile/designfile.tsv' --comparefile 'path/to/comparefile/compare.txt' -profile <docker/test/conda> -resume

    Mandatory arguments:
      --designfile                  .tsv format(table splited by tabs): Sample_ID, input_R1, input_R2, ip_R1, ip_R2, Group_ID
      --comparefile                 .txt format: control_vs_treated;
                                    '--comparefile false' for projects without differential analysis
      -profile                      Configuration files can contain the definition of one or more profiles.
                                    Available: conda, docker, test
    
    References                      If not specified in the configuration file or you wish to overwrite any of the references.
      --fasta                       Path to Fasta reference
      --gtf                         Path to GTF reference

    Index:
      --rRNA_fasta                  Path to HISAT2 rRNA fasta reference
      --tophat2_index               Path to TopHat2 index, eg. 'path/to/Tophat2Index/*'
      --hisat2_index                Path to HISAT2 index, eg. 'path/to/HISAT2Index/*'
      --bwa_index                   Path to BWA index, eg. 'path/to/BWAIndex/*'
      --star_index                  Path to STAR index, eg. 'path/to/STARIndex/*'

    Options:
      --skip_qc                     Skip all Quality Control steps                        
      --skip_peakCalling            Skip all Peak Calling steps
      --skip_diffpeakCalling        Skip all Differential methylation analysis
      --skip_fastqc                 Skip FastQC
      --skip_rseqc                  Skip RSeQC
      --skip_genebody_coverage      Skip calculating genebody coverage  
      --skip_edger                  Skip the EdgeR process of differential expression analysis steps
      --skip_deseq2                 Skip the DESeq2 process of differential expression analysis steps
      --skip_metpeak                Skip the MeTPeak process of Peak Calling steps
      --skip_macs2                  Skip the MACS2 process of Peak Calling steps
      --skip_matk                   Skip the MATK process of Peak Calling steps
      --skip_QNB                    Skip the QNB process of Differential methylation analysis
      --skip_diffmatk               Skip the MATK process of Differential methylation analysis
    
    Main parameters of analysis mode:
      --stranded                    "yes" OR "no" OR "reverse"
      --mapq_cutoff                 [0-255] for filtering reads of different mapping quality
      --aligners                    "star" OR "bwa" OR "tophat2" OR "hisat2" OR "none","none" for bam files
      --expression_analysis_mode    "DESeq2" OR "edgeR" OR "none"
      --peakCalling_mode            "group" OR "independence" for MATK and MeTPeak
      --peakMerged_mode             "rank" OR "macs2" OR "MATK" OR "metpeak" OR "mspc"
      --methylation_analysis_mode   "MATK" OR "QNB" OR "Wilcox-test" OR "edgeR" OR "DESeq2"

    Other options:
      --outdir                      The output directory where the results will be saved, defalut = $baseDir/results
      --email                       Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits
      --help                        To show the help message
      -name                         Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.

    AWSBatch options:
      --awsqueue                    The AWSBatch JobQueue that needs to be set when running on AWSBatch
      --awsregion                   The AWS Region for your AWS Batch job to run on
    """.stripIndent()
}

// Show help message
if (params.help) {
    helpMessage()
    exit 0
}

// Stage config files
ch_multiqc_config = file(params.multiqc_config, checkIfExists: true)
ch_output_docs = file("$baseDir/docs/output.md", checkIfExists: true)

/*
 * SET UP CONFIGURATION VARIABLES
 */

// Configurable variables
params.name = false
params.project = false
params.genome = false
params.call = false
params.email = false
params.plaintext_email = false
params.seqCenter = false
params.help = false

// Validate inputs
// Check if genome exists in the config file
if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
    exit 1, "The provided genome '${params.genome}' is not available in the iGenomes file. Currently the available genomes are ${params.genomes.keySet().join(", ")}"
}

params.star_index = params.genome ? params.genomes[ params.genome ].star ?: false : false
params.fasta = params.genome ? params.genomes[ params.genome ].fasta ?: false : false
params.gtf = params.genome ? params.genomes[ params.genome ].gtf ?: false : false
params.bed12 = params.genome ? params.genomes[ params.genome ].bed12 ?: false : false
params.hisat2_index = params.genome ? params.genomes[ params.genome ].hisat2 ?: false : false

// Has the run name been specified by the user?
// this has the bonus effect of catching both -name and --name
custom_runName = params.name
if (!(workflow.runName ==~ /[a-z]+_[a-z]+/)) {
  custom_runName = workflow.runName
}

if ( workflow.profile == 'awsbatch') {
  // AWSBatch sanity checking
  if (!params.awsqueue || !params.awsregion) exit 1, "Specify correct --awsqueue and --awsregion parameters on AWSBatch!"
  // Check outdir paths to be S3 buckets if running on AWSBatch
  // related: https://github.com/nextflow-io/nextflow/issues/813
  if (!params.outdir.startsWith('s3:')) exit 1, "Outdir not on S3 - specify S3 Bucket to run on AWSBatch!"
  // Prevent trace files to be stored on S3 since S3 does not support rolling files.
  if (workflow.tracedir.startsWith('s3:')) exit 1, "Specify a local tracedir or run without trace! S3 cannot be used for tracefiles."
}

//check fasta files
if ( params.fasta ){
    fasta = file(params.fasta, checkIfExists: true)
    if( !fasta.exists() ) exit 1, LikeletUtils.print_red("Fasta file not found: ${params.fasta}")
}else {
    exit 1, LikeletUtils.print_red("No reference genome specified!")
}

//check gtf files
if( params.gtf ){
    gtf = file ( params.gtf, checkIfExists: true)
    if( !gtf.exists() ) exit 1, LikeletUtils.print_red("gtf not found: ${params.gtf}")
} else {
    exit 1, LikeletUtils.print_red("No GTF annotation specified!")
}

//check rRNA fasta files
if ( params.rRNA_fasta ){
    rRNA_fasta = file(params.rRNA_fasta, checkIfExists: true)
    if( !rRNA_fasta.exists() ) exit 1, LikeletUtils.print_red("rRNA fasta file not found: ${params.rRNA_fasta}")
}else{
    rRNA_fasta = Channel.from("")
}

//check comparefile
if( params.comparefile == "two_group" ){
    comparefile = false
    compareLines = Channel.from("two_group")
}else if( params.comparefile){
    comparefile = file(params.comparefile)
    if( !comparefile.exists() ) exit 1, print_red("Compare file not found: ${params.comparefile}")
    compareLines = Channel.from(comparefile.readLines())
}else {
    comparefile = false
    compareLines = Channel.from("")
}

compareLines.into{
    compareLines_for_DESeq2; compareLines_for_edgeR; compareLines_for_plot;
    compareLines_for_diffm6A; compareLines_for_arranged_result
}

// Validate the params of skipping Aligners Tools Setting
if( params.aligners == "none" || params.aligners == "star" || params.aligners == "hisat2" || params.aligners == "tophat2" || params.aligners == "bwa" ){
    aligner = params.aligners
}else{
    exit 1, LikeletUtils.print_red("Invalid aligner option: ${params.aligner}. Valid options: 'star', 'hisat2', 'tophat2', 'bwa', 'none'")
}
if( params.expression_analysis_mode == "edgeR" ){
    params.skip_edger = false
    params.skip_deseq2 = true
    params.skip_expression = false
}else if( params.expression_analysis_mode == "DESeq2" ){
    params.skip_edger = true
    params.skip_deseq2 = false
    params.skip_expression = false
}else if( params.expression_analysis_mode == "none" ){
    params.skip_edger = true
    params.skip_deseq2 = true
    params.skip_expression = true
}else{
    exit 1, LikeletUtils.print_red("Invalid expression_analysis_mode option: ${params.expression_analysis_mode}. Valid options: 'edgeR', 'DESeq2', 'none'")
}


/*
 * Create a channel for input read files
*/
if ( params.readPaths ){
    if (aligner == 'none') {
        Channel
            .fromPath( params.readPaths )
            .ifEmpty { exit 1, "params.readPaths was empty - no input files supplied" }
            .into{ raw_fastq; raw_bam }
    } else if ( params.single_end ) {
        Channel
            .from( params.readPaths )
            .map{ row -> [ row[0], [ file(row[1][0], checkIfExists: true) ] ] }
            .ifEmpty { exit 1, "params.readPaths was empty - no input files supplied" }
            .into{ raw_fastq; raw_bam }
    } else {
        Channel
            .from( params.readPaths )
            .map{ row -> [ row[0], [ file(row[1][0], checkIfExists: true), file(row[1][1], checkIfExists: true) ] ] }
            .ifEmpty { exit 1, "params.readPaths was empty - no input files supplied" }
            .into{ raw_fastq; raw_bam }
    }
}else if( params.designfile ) {
    designfile = file(params.designfile, checkIfExists: true)
    if( !designfile.exists() ) exit 1, LikeletUtils.print_red("Design file not found: ${params.designfile}")
    LikeletUtils.extractData(designfile).into{ input_data_fastq; input_data_bam; designinfo}
    input_data_fastq.filter{ it[6] == "fastq" }.map{ it.subList(0,6) }.set{ raw_fastq } //filter filetype & delete filetype
    input_data_bam.filter{ it[6] == "bam" }.map{ it.subList(0,6) }.set{ raw_bam }
    designinfo.map{ it.getAt([0,5]) }.set{ format_design }
}else{
    exit 1, LikeletUtils.print_red("No Design file specified!")
}

/*
                         showing the process and files
========================================================================================
*/
// Header log info	
println LikeletUtils.sysucc_ascii()
println LikeletUtils.print_purple("============You are running MeRIPseqPipe with the following parameters===============")
println LikeletUtils.print_purple("Checking parameters ...")
println LikeletUtils.print_yellow("===================================Pipeline summary=============================")
println (LikeletUtils.print_yellow("Max Resources                  : ") + LikeletUtils.print_green("$params.max_memory memory, $params.max_cpus cpus, $params.max_time time per job"))
if (workflow.containerEngine){
    println (LikeletUtils.print_yellow("Container                       : ") + LikeletUtils.print_green("$workflow.containerEngine - $workflow.container"))
}
println (LikeletUtils.print_yellow("Output dir                     : ") + LikeletUtils.print_green(params.outdir))
println (LikeletUtils.print_yellow("Launch dir                     : ") + LikeletUtils.print_green(workflow.launchDir))
println (LikeletUtils.print_yellow("Working dir                    : ") + LikeletUtils.print_green(workflow.workDir))
println (LikeletUtils.print_yellow("Script dir                     : ") + LikeletUtils.print_green(workflow.projectDir))
println (LikeletUtils.print_yellow("User                           : ") + LikeletUtils.print_green(workflow.userName))
if (workflow.profile == 'awsbatch') {
    println (LikeletUtils.print_yellow("AWS Region                      : ") + LikeletUtils.print_green(params.awsregion))
    println (LikeletUtils.print_yellow("AWS Queue                       : ") + LikeletUtils.print_green(params.awsqueue))
}
println (LikeletUtils.print_yellow("Config Profile                 : ") + LikeletUtils.print_green(workflow.profile))
if (params.config_profile_description){
    println (LikeletUtils.print_yellow("Config Description              : ") + LikeletUtils.print_green(params.config_profile_description))
}
if (params.config_profile_contact){
    println (LikeletUtils.print_yellow("Config Contact                  : ") + LikeletUtils.print_green(params.config_profile_contact))
}
if (params.config_profile_url){
    println (LikeletUtils.print_yellow("Config URL                      : ") + LikeletUtils.print_green(params.config_profile_url))
}
if (params.email || params.email_on_fail) {
    println (LikeletUtils.print_yellow("E-mail Address                  : ") + LikeletUtils.print_green(params.email))
    println (LikeletUtils.print_yellow("E-mail on failure               : ") + LikeletUtils.print_green(params.email_on_fail))
    println (LikeletUtils.print_yellow("MultiQC maxsize                 : ") + LikeletUtils.print_green(params.maxMultiqcEmailFileSize))
}
println LikeletUtils.print_yellow("=====================================Reads types================================")
println (LikeletUtils.print_yellow("SingleEnd                      : ") + LikeletUtils.print_green(params.single_end ? 'Single-End' : 'Paired-End'))
println (LikeletUtils.print_yellow("Stranded                       : ") + LikeletUtils.print_green(params.stranded))

println LikeletUtils.print_yellow("====================================Mode selected==============================")
println (LikeletUtils.print_yellow("aligners                       : ") + LikeletUtils.print_green(params.aligners))
println (LikeletUtils.print_yellow("peakCalling_mode               : ") + LikeletUtils.print_green(params.peakCalling_mode))
println (LikeletUtils.print_yellow("peakMerged_mode                : ") + LikeletUtils.print_green(params.peakMerged_mode))
println (LikeletUtils.print_yellow("expression_analysis_mode       : ") + LikeletUtils.print_green(params.expression_analysis_mode))
println (LikeletUtils.print_yellow("methylation_analysis_mode      : ") + LikeletUtils.print_green(params.methylation_analysis_mode))

println LikeletUtils.print_yellow("==================================Input files selected==========================")
println (LikeletUtils.print_yellow("Reads Path                     : ") + LikeletUtils.print_green(params.readPaths ? "github" : params.input))
println (LikeletUtils.print_yellow("fasta file                     : ") + LikeletUtils.print_green(params.fasta))
println (LikeletUtils.print_yellow("Gtf file                       : ") + LikeletUtils.print_green(params.gtf))
println (LikeletUtils.print_yellow("Design file                    : ") + LikeletUtils.print_green(params.designfile))
println (LikeletUtils.print_yellow("Compare file                   : ") + LikeletUtils.print_green(params.comparefile))

println LikeletUtils.print_yellow("==================================Skip model selected==========================")
println (LikeletUtils.print_yellow("Skip samtools sort             : ") + LikeletUtils.print_green(params.skip_sort))
println (LikeletUtils.print_yellow("Skip expression analysis       : ") + LikeletUtils.print_green(params.skip_expression))
println (LikeletUtils.print_yellow("Skip peakCalling               : ") + LikeletUtils.print_green(params.skip_peakCalling))
println (LikeletUtils.print_yellow("Skip diffpeakCalling           : ") + LikeletUtils.print_green(params.skip_diffpeakCalling))
println (LikeletUtils.print_yellow("Skip annotation                : ") + LikeletUtils.print_green(params.skip_annotation))
println (LikeletUtils.print_yellow("Skip qc                        : ") + LikeletUtils.print_green(params.skip_qc))

// log.info summary.collect { k,v -> "${k.padRight(18)}: $v" }.join("\n")
// // Check the hostnames against configured profiles
// checkHostname()
/*