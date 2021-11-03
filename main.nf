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
 * Bao Xiaoqiong <baoxq@sysucc.org.cn>
 * Zhu Kaiyu <zhuky5@mail2.sysu.edu.cn>:
 * Qi Zhao <zhaoqi@sysucc.org.cn>;
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

    nextflow run path/to/MeRIPseqPipe --designfile 'path/to/designfile.tsv' --comparefile 'path/to/compare.txt' -profile <docker/test/conda> -resume
    
    OR, 

    nextflow path/to/MeRIPseqPipe/main.nf -c path/to/MeRIPseqPipe/nextflow.config --designfile 'path/to/designfile/designfile.tsv' --comparefile 'path/to/comparefile/compare.txt' -profile <docker/test/conda> -resume

    Mandatory arguments:
      --designfile                  .tsv format(table splited by tabs): Sample_ID, input_R1, input_R2, ip_R1, ip_R2, Group_ID
      --comparefile                 .txt format: control_vs_treated;
                                    "false" for projects without differential analysis and "two_groups" for only two groups in the designfile
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

    Main parameters of analysis mode:
      --stranded                    "yes" OR "no" OR "reverse"
      --mapq_cutoff                 [0-255] for filtering reads of different mapping quality
      --featurecount_minMQS         Integer giving the minimum mapping quality score a read must satisfy in order to be counted.
      --motiflength                 length for HOMER motif searching
      --aligners                    "star" OR "bwa" OR "tophat2" OR "hisat2" OR "none","none" for bam files
      --expression_analysis_mode    "DESeq2" OR "edgeR" OR "none"
      --peakCalling_mode            "group" OR "independence" for MATK and MeTPeak
      --peakMerged_mode             "rank" OR "macs2" OR "MATK" OR "metpeak" OR "mspc"
      --methylation_analysis_mode   "MATK" OR "QNB" OR "Wilcox-test" OR "edgeR" OR "DESeq2"

    Process skipping options:
      --skip_fastp                  Skip fastp
      --skip_fastqc                 Skip FastQC
      --skip_rseqc                  Skip RSeQC
      --skip_filterrRNA             Skip the process of rRNA cleaning before mapping
      --skip_qc                     Skip all Quality Control steps
      --skip_edger                  Skip the EdgeR process of differential expression analysis steps
      --skip_deseq2                 Skip the DESeq2 process of differential expression analysis steps
      --skip_sort                   Skip the process of sorting BAM files
      --skip_expression             Skip all Differential expression analysis
      --skip_peakCalling            Skip all Peak Calling steps
      --skip_diffpeakCalling        Skip all Differential methylation analysis
      --skip_metpeak                Skip the MeTPeak process of Peak Calling steps
      --skip_macs2                  Skip the MACS2 process of Peak Calling steps
      --skip_matk                   Skip the MATK process of Peak Calling steps
      --skip_meyer                  Skip the Meyer process of Peak Calling steps
      --skip_QNB                    Skip the QNB process of Differential methylation analysis
      --skip_diffmatk               Skip the MATK process of Differential methylation analysis
      --skip_annotation             Skip the process of peak annotation
      --skip_motif                  Skip the process of motif searching
      --skip_m6Aprediction          Skip the process of m6A site prediction
      --skip_createbedgraph         Skip the process of making bedgraph format files and creatIGVjs
    
    Other options:
      --outdir                      The output directory where the results will be saved, defalut = $baseDir/results
      --help                        To show the help message
      -name                         Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.

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
if( params.comparefile == "two_groups" ){
    comparefile = false
    compareLines = Channel.from("two_groups")
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
println (LikeletUtils.print_yellow("gzip                           : ") + LikeletUtils.print_green(params.gzip))

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

========================================================================================
                             check or build the index
========================================================================================
*/ 
/*
 * PREPROCESSING - Build BED12 file
 * NEED gtf.file
 */
 process CheckDesignCompare{
    input:
    val design_info from format_design.collect()
    file comparefile

    output:
    file (formatted_design) into formatted_designfile

    script:
    formatted_design = "formatted_designfile.txt"
    formatted_design_info = ""
    for(int i = 0; i < design_info.size(); i+=2 ) {
        sample = design_info[i] + ".input," + design_info[i] + ".ip"
        formatted_design_info += design_info[i] + "," + sample + "," + design_info[i+1] + "\n"
    }
    """
    echo "Sample_ID,input_FileName,ip_FileName,Group" > $formatted_design
    echo "$formatted_design_info" |awk NF |sort | uniq >> $formatted_design
    # Check the consistency of designfile and comparefile
    if [ "$comparefile" != "input.1" ]; then 
        ## get groups' name in comparefile
        cat $comparefile | dos2unix | awk -F "_vs_" '{print \$1"\\n"\$2}' | sort | uniq > tmp.compare.group
        ## get groups' name in designfile
        awk -F, 'NR>1{print \$4}' $formatted_design | sort | uniq > tmp.design.group
        intersection_num=\$(join tmp.compare.group tmp.design.group | wc -l)
        if [[ \$intersection_num  != \$(cat tmp.compare.group| wc -l) ]] ;then 
            echo "The groups' name of comparefile and designfile are inconsistent."
            echo "Please check your comparefile: "$comparefile
            echo "The groups' name of comparefile in designfile: "\$(join tmp.compare.group tmp.design.group)
            exit 1
        fi
        rm tmp.compare.group tmp.design.group
    fi
    """
}


process makeBED12 {
    label 'build_index'
    tag "gtf2bed12"
    publishDir path: { params.saveReference ? "${params.outdir}/Genome/reference_genome" : params.outdir },
                saveAs: { params.saveReference ? it : null }, mode: 'copy'

    when:
    !params.skip_qc && !params.skip_rseqc || !params.skip_motif

    input:
    file gtf
    
    output:
    file "${gtf.baseName}.bed" into bed12file

    shell:      
    '''
    gtf_file=!{gtf}
    gtfToGenePred -genePredExt -geneNameAsName2 $gtf_file ${gtf_file/.gtf/.tmp}
    genePredToBed ${gtf_file/.gtf/.tmp} ${gtf_file/.gtf/.bed}
    '''
}

process makechromesize {
    label 'build_index'
    tag "gtf2bed12"
    publishDir path: { params.saveReference ? "${params.outdir}/Genome/reference_genome" : params.outdir },
                saveAs: { params.saveReference ? it : null }, mode: 'copy'

    when:
    true

    input:
    file fasta
    
    output:
    file "chromsizes.file" into chromsizesfile

    shell:      
    """
    samtools faidx ${fasta}
    cut -f1,2 ${fasta}.fai > chromsizes.file
    """
}

/*
 * PREPROCESSING - Build TOPHAT2 index
 * NEED genome.fa
 */
if( params.tophat2_index && aligner == "tophat2" ){
    tophat2_index = Channel
        .fromPath( params.tophat2_index )
        .ifEmpty { exit 1, "Tophat2 index not found: ${params.tophat2_index}" }
}else if( params.fasta ){
    process MakeTophat2Index {
        label 'build_index'
        tag "tophat2_index"
        publishDir path: { params.saveReference ? "${params.outdir}/Genome/": params.outdir },
                   saveAs: { params.saveReference ? it : null }, mode: 'copy'
        input:
        file fasta

        output:
        file "Tophat2Index/*" into tophat2_index

        when:
        aligner == "tophat2"

        script:
        tophat2_index = "Tophat2Index/" + fasta.baseName.toString()
        """
        mkdir Tophat2Index
        ln $fasta Tophat2Index
        bowtie2-build -f $fasta $tophat2_index
        """
    }
}else {
    exit 1, println LikeletUtils.print_red("There is no Tophat2 Index")
}

/*
 * PREPROCESSING - Build HISAT2 index
 * NEED genome.fa genes.gtf snp.txt/vcf
 */
if( params.hisat2_index && aligner == "hisat2" ){
    hisat2_index = Channel
        .fromPath(params.hisat2_index)
        .ifEmpty { exit 1, "hisat2 index not found: ${params.hisat2_index}" }
}else if( params.fasta){
    process MakeHisat2Index {
        label 'build_index'
        tag "hisat2_index"
        publishDir path: { params.saveReference ? "${params.outdir}/Genome/ " : params.outdir },
                   saveAs: { params.saveReference ? it : null }, mode: 'copy'        
        input:
        file fasta
        file gtf

        output:
        file "Hisat2Index/*" into hisat2_index

        when:
        aligner == "hisat2"
        
        script:
        """
        mkdir Hisat2Index
        hisat2_extract_exons.py $gtf > Hisat2Index/${fasta.baseName}.exon
        hisat2_extract_splice_sites.py $gtf > Hisat2Index/${fasta.baseName}.ss
        hisat2-build -p ${task.cpus} -f $fasta --exon Hisat2Index/${fasta.baseName}.exon --ss Hisat2Index/${fasta.baseName}.ss Hisat2Index/${fasta.baseName}
        """
    }
}else {
    exit 1, println LikeletUtils.print_red("There is no Hisat2 Index")
}

/*
 * PREPROCESSING - Build BWA index
 * NEED genome.fa
 */
if( params.bwa_index && aligner == "bwa" ){
    bwa_index = Channel
        .fromPath( params.bwa_index )
        .ifEmpty { exit 1, "bwa index not found: ${params.bwa_index}" }
}else if(params.fasta ){
    process MakeBWAIndex {
        label 'build_index'
        tag "bwa_index"
        publishDir path: { params.saveReference ? "${params.outdir}/Genome/" : params.outdir },
                   saveAs: { params.saveReference ? it : null }, mode: 'copy'

        input:
        file fasta

        output:
        file "BWAIndex/*" into bwa_index

        when:
        aligner == "bwa"
     
        script:
        """
        mkdir BWAIndex
        cd BWAIndex/
        bwa index -p ${fasta.baseName} -a bwtsw ../$fasta
        cd ../
        """
    }
}else {
    exit 1, println LikeletUtils.print_red("There is no BWA Index")
}

/*
 * PREPROCESSING - Build STAR index
 * NEED genome.fa genes.gtf
 */
if( params.star_index && aligner == "star"){
    star_index = Channel
        .fromPath(params.star_index)
        .ifEmpty { exit 1, "STAR index not found: ${params.star_index}" }
}else if( params.fasta ){
    process MakeStarIndex {
        label 'build_index'
        tag "star_index"
        publishDir path: { params.saveReference ? "${params.outdir}/Genome/" : params.outdir },
                   saveAs: { params.saveReference ? it : null }, mode: 'copy'
        input:
        file fasta
        file gtf

        output:
        file "StarIndex" into star_index

        when:
        aligner == "star"

        script:
        readLength = 50
        overhang = readLength - 1
        """
        mkdir StarIndex
        STAR --runThreadN ${task.cpus} \
        --runMode genomeGenerate \
        --genomeDir StarIndex \
        --genomeFastaFiles $fasta \
        --sjdbGTFfile $gtf \
        --sjdbOverhang $overhang \
	--limitGenomeGenerateRAM 36000000000
        """
    }
}else {
   exit 1, println LikeletUtils.print_red("There is no STAR Index")
}

/*
 * PREPROCESSING - Build rRNA index
 * NEED rRNA.fa 
 */
process MakerRNAindex {
    label 'build_index'
    tag "rRNA_index"
    publishDir path: { params.saveReference ? "${params.outdir}/Genome/" : params.outdir },
                saveAs: { params.saveReference ? it : null }, mode: 'copy'
    input:
    file rRNA_fasta from rRNA_fasta

    output:
    file "rRNAindex/*" into rRNA_index

    when:
    params.rRNA_fasta && !params.skip_filterrRNA

    script:
    """
    mkdir rRNAindex
    hisat2-build -p ${task.cpus} -f $rRNA_fasta rRNAindex/${rRNA_fasta.baseName}
    """
}

/*
========================================================================================
                                Step 1. QC------Fastp & FastQC
========================================================================================
*/ 
process Fastp{
    label 'fastp'
    tag "$sample_name"
    //errorStrategy 'ignore'
    publishDir path: { params.skip_fastp ? params.outdir : "${params.outdir}/QC/fastp" },
             saveAs: { params.skip_fastp ? null : it }, mode: 'link'
        
    input:
    set val(sample_id), file(reads), val(reads_single_end), val(gzip), val(input), val(group) from raw_fastq

    output:
    set val(sample_name), file("*_aligners.fastq*"), val(reads_single_end), val(sample_id), val(gzip), val(input), val(group) into fastqc_reads, fastp_reads
    file "*" into fastp_results

    when:
    aligner != "none"

    shell:
    skip_fastp = params.skip_fastp
    if ( reads_single_end ){
        filename = reads.toString() - ~/(\.fq)?(\.fastq)?(\.gz)?$/
        sample_name = filename
        add_aligners = sample_name + "_aligners.fastq" + (gzip ? ".gz" : "")
        """
        if [ $skip_fastp == "false" ]; then
            fastp -i ${reads} -o ${add_aligners} -j ${sample_name}_fastp.json -h ${sample_name}_fastp.html -w ${task.cpus}
        else
            mv ${reads} ${add_aligners}
        fi
        """
    } else {
        filename = reads[0].toString() - ~/(_R[0-9])?(_[0-9])?(\.fq)?(\.fastq)?(\.gz)?$/
        sample_name = filename
        add_aligners_1 = sample_name + "_1_aligners.fastq" + (gzip ? ".gz" : "")
        add_aligners_2 = sample_name + "_2_aligners.fastq" + (gzip ? ".gz" : "")
        """
        if [ $skip_fastp == "false" ]; then  
            fastp -i ${reads[0]} -o ${add_aligners_1} -I ${reads[1]} -O ${add_aligners_2} -j ${sample_name}_fastp.json -h ${sample_name}_fastp.html -w ${task.cpus}
        else
            mv ${reads[0]} ${add_aligners_1}
            mv ${reads[1]} ${add_aligners_2}
        fi
        """
    } 
}

process Fastqc{
    tag "$sample_name"
    publishDir path: { params.skip_fastqc ? params.outdir : "${params.outdir}/QC" },
             saveAs: { params.skip_fastqc ? null : it }, mode: 'link'

    input:
    set val(sample_name), file(reads), val(reads_single_end), val(sample_id), val(gzip), val(input), val(group) from fastqc_reads

    output:
    file "fastqc/*" into fastqc_results

    when:
    aligner != "none" && !params.skip_fastqc

    shell:
    skip_fastqc = params.skip_fastqc
    if ( reads_single_end){
        """
        mkdir fastqc
        fastqc -o fastqc --noextract ${reads}
        """       
    } else {
        """
        mkdir fastqc   
        fastqc -o fastqc --noextract ${reads[0]}
        fastqc -o fastqc --noextract ${reads[1]}
        """      
    }
}

/*
========================================================================================
                            Step 2. Reads Mapping
========================================================================================
*/

if(params.rRNA_fasta && !params.skip_filterrRNA){
    fastp_reads.set{rRNA_reads}

    process FilterrRNA {
    label 'aligners'
    tag "$sample_name"
    publishDir "${params.outdir}/alignment/rRNA_dup", mode: 'link', overwrite: true
    
    input:
    set val(sample_name), file(reads), val(reads_single_end), val(sample_id), val(gzip), val(input), val(group) from rRNA_reads
    file index from rRNA_index.collect()

    output:
    set val(sample_name), file("*.fastq.gz"), val(reads_single_end), val(sample_id), val(gzip), val(input), val(group) into tophat2_reads, hisat2_reads, bwa_reads, star_reads
    file "*_summary.txt" into rRNA_log

    when:
    params.rRNA_fasta && !params.skip_filterrRNA

    script:
    gzip = true
    index_base = index[0].toString() - ~/(\.exon)?(\.\d)?(\.fa)?(\.gtf)?(\.ht2)?$/
    if (reads_single_end) {
        """
        hisat2 --summary-file ${sample_name}_rRNA_summary.txt \
            --no-spliced-alignment --no-softclip --norc --no-unal \
            -p ${task.cpus} --dta --un-gz ${sample_id}.fastq.gz \
            -x $index_base \
            -U $reads | \
            samtools view -@ ${task.cpus} -Shub - | \
            samtools sort -@ ${task.cpus} -o ${sample_name}_rRNA_sort.bam -
        """
    } else {
        """
        hisat2 --summary-file ${sample_name}_rRNA_summary.txt \
            --no-spliced-alignment --no-softclip --norc --no-unal \
            -p ${task.cpus} --dta --un-conc-gz ${sample_name}_fastq.gz \
            -x $index_base \
            -1 ${reads[0]} -2 ${reads[1]} | \
            samtools view -@ ${task.cpus} -Shub - | \
            samtools sort -@ ${task.cpus} -o ${sample_name}_rRNA_sort.bam -
        mv ${sample_name}_fastq.1.gz ${sample_name}_1.fastq.gz
        mv ${sample_name}_fastq.2.gz ${sample_name}_2.fastq.gz
        """
    }
    }
}else{
    fastp_reads.into{tophat2_reads; hisat2_reads; bwa_reads; star_reads}
}

process Tophat2Align {
    label 'aligners'
    tag "$sample_name"
    publishDir "${params.outdir}/alignment/tophat2", mode: 'link', overwrite: true

    input:
    set val(sample_name), file(reads), val(reads_single_end), val(sample_id), val(gzip), val(input), val(group) from tophat2_reads
    file index from tophat2_index.collect()
    file gtf

    output:
    set val(sample_id), file("*_tophat2.bam"), val(reads_single_end), val(gzip), val(input), val(group) into tophat2_bam
    file "*_log.txt" into tophat2_log
    
    when:
    aligner == "tophat2"

    script:
    index_base = index[0].toString() - ~/(\.rev)?(\.\d)?(\.fa)?(\.bt2)?$/
    strand_info = params.stranded == "no" ? "fr-unstranded" : params.stranded == "reverse" ? "fr-secondstrand" : "fr-firststrand"
    if (reads_single_end) {
        """
        tophat  -p ${task.cpus} \
                -G $gtf \
                -o $sample_name \
                --no-novel-juncs \
                --library-type $strand_info \
                $index_base \
                $reads > ${sample_name}_log.txt
        mv $sample_name/accepted_hits.bam ${sample_name}_tophat2.bam
        """
    } else {
        """
        tophat -p ${task.cpus} \
                -G $gtf \
                -o $sample_name \
                --no-novel-juncs \
                --library-type $strand_info \
                $index_base \
                ${reads[0]} ${reads[1]} > ${sample_name}_log.txt
        mv $sample_name/accepted_hits.bam ${sample_name}_tophat2.bam
        """
    }
}

process Hisat2Align {
    label 'aligners'
    tag "$sample_name"
    publishDir "${params.outdir}/alignment/hisat2", mode: 'link', overwrite: true

    input:
    set val(sample_name), file(reads), val(reads_single_end), val(sample_id), val(gzip), val(input), val(group) from hisat2_reads
    file index from hisat2_index.collect()

    output:
    set val(sample_id), file("*_hisat2.bam"), val(reads_single_end), val(gzip), val(input), val(group) into hisat2_bam
    file "*_summary.txt" into hisat2_log

    when:
    aligner == "hisat2"

    script:
    index_base = index[0].toString() - ~/(\.exon)?(\.\d)?(\.fa)?(\.gtf)?(\.ht2)?$/
    if (reads_single_end) {
        """
        hisat2  --summary-file ${sample_name}_hisat2_summary.txt\
                -p ${task.cpus} --dta \
                -x $index_base \
                -U $reads | \
                samtools view -@ ${task.cpus} -hbS - > ${sample_name}_hisat2.bam 
        """
    } else {
        """
        hisat2  --summary-file ${sample_name}_hisat2_summary.txt \
                -p ${task.cpus} --dta \
                -x $index_base \
                -1 ${reads[0]} -2 ${reads[1]} | \
                samtools view -@ ${task.cpus} -hbS - > ${sample_name}_hisat2.bam
        """
    }
}

process BWAAlign{
    label 'aligners'
    tag "$sample_name"
    publishDir "${params.outdir}/alignment/bwa", mode: 'link', overwrite: true
    
    input:
    set val(sample_name), file(reads), val(reads_single_end), val(sample_id), val(gzip), val(input), val(group) from bwa_reads
    file index from bwa_index.collect()

    output:
    set val(sample_id), file("*_bwa.bam"), val(reads_single_end), val(gzip), val(input), val(group) into bwa_bam
    file "*" into bwa_result

    when:
    aligner == "bwa"

    script:
    index_base = index[0].toString() - ~/(\.pac)?(\.bwt)?(\.ann)?(\.amb)?(\.sa)?(\.fa)?$/
    if (reads_single_end) {
        """
        bwa aln -t ${task.cpus} \
                -f ${reads.baseName}.sai \
                $index_base \
                $reads
        bwa samse \
                $index_base \
                ${reads.baseName}.sai \
                $reads  | \
            samtools view -@ ${task.cpus} -h -bS - > ${sample_name}_bwa.bam
        """
    } else {
        """
        bwa aln -t ${task.cpus} \
                -f ${reads[0].baseName}.sai \
                $index_base \
                ${reads[0]}
        bwa aln -t ${task.cpus} \
                -f ${reads[1].baseName}.sai \
                $index_base \
                ${reads[1]}
        bwa sampe \
                $index_base \
                ${reads[0].baseName}.sai ${reads[1].baseName}.sai \
                ${reads[0]} ${reads[1]} | \
            samtools view -@ ${task.cpus} -h -bS - > ${sample_name}_bwa.bam
        """
    }
}

process StarAlign {
    label 'aligners'
    tag "$sample_name"
    publishDir "${params.outdir}/alignment/star", mode: 'link', overwrite: true
    
    input:
    set val(sample_name), file(reads), val(reads_single_end), val(sample_id), val(gzip), val(input), val(group) from star_reads
    file star_index from star_index.collect()

    output:
    set val(sample_id), file("*_star.bam"), val(reads_single_end), val(gzip), val(input), val(group) into star_bam
    file "*.final.out" into star_log

    when:
    aligner == "star"

    script:
    gzip_cmd = gzip ? "--readFilesCommand zcat" : ""
    if (reads_single_end) {
        """
        STAR --runThreadN ${task.cpus} $gzip_cmd \
            --twopassMode Basic \
            --genomeDir $star_index \
            --readFilesIn $reads  \
            --outSAMtype BAM Unsorted \
            --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 \
            --outFilterIntronMotifs RemoveNoncanonical \
            --outFilterMultimapNmax 20 \
            --alignIntronMin 20 \
            --alignIntronMax 100000 \
            --alignMatesGapMax 1000000 \
            --outFileNamePrefix ${sample_name}  > ${sample_name}_log.txt
        mv ${sample_name}Aligned.out.bam ${sample_name}_star.bam
        """
    } else {
        """
        STAR --runThreadN ${task.cpus} $gzip_cmd \
            --twopassMode Basic \
            --genomeDir $star_index \
            --readFilesIn ${reads[0]} ${reads[1]}  \
            --outSAMtype BAM Unsorted \
            --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 \
            --outFilterIntronMotifs RemoveNoncanonical \
            --outFilterMultimapNmax 20 \
            --alignIntronMin 20 \
            --alignIntronMax 1000000 \
            --alignMatesGapMax 1000000 \
            --outFileNamePrefix ${sample_name} > ${sample_name}_log.txt
        mv ${sample_name}Aligned.out.bam ${sample_name}_star.bam
        """
    }
}


/*
========================================================================================
                        Step 3 Sort BAM file AND QC
========================================================================================
*/ 
Channel
    .from()
    .concat(tophat2_bam, hisat2_bam, bwa_bam, star_bam, raw_bam)
    .set{ merge_bam_file }
/*
 * STEP 3-1 - Sort BAM file
*/

process SortRename {
    label 'sort'
    tag "$sample_name"
    publishDir "${params.outdir}/alignment/samtoolsSort/", mode: 'link', overwrite: true
    
    input:
    set val(sample_id), file(bam_file), val(reads_single_end), val(gzip), val(input), val(group) from merge_bam_file

    output:
    set val(group), val(sample_id), file("*.bam"), file("*.bai") into sorted_bam
    file "*.{bam,bai}" into rseqc_bam, bedgraph_bam, feacount_bam, cuffbam, peakquan_bam, diffpeak_bam, sng_bam
    file "*.bam" into bam_results

    script:
    sample_name = sample_id + (input ? ".input_" : ".ip_") + group
    output = sample_name + ".bam"
    mapq_cutoff = (params.mapq_cutoff).toInteger() 
    if (!params.skip_sort){
        """
        if [ "$mapq_cutoff" -gt "0" ]; then
            samtools view -hbq $mapq_cutoff $bam_file | samtools sort -@ ${task.cpus} -O BAM -o $output -
        else
            samtools sort -@ ${task.cpus} -O BAM -o $output $bam_file
        fi
        samtools index -@ ${task.cpus} $output
        """
    } else {
        """
        if [ "$mapq_cutoff" -gt "0" ]; then
            samtools view -hbq $mapq_cutoff $bam_file > $output
        else
            mv $bam_file $output
        fi
        samtools index -@ ${task.cpus} $output
        """
    }
}


/*
 * STEP 3-2 - RSeQC analysis
*/
process RSeQC {
    tag "$output"
    publishDir "${params.outdir}/QC/rseqc" , mode: 'link', overwrite: true,
        saveAs: {filename ->
                 if (filename.indexOf("bam_stat.txt") > 0)                      "bam_stat/$filename"
            else if (filename.indexOf("infer_experiment.txt") > 0)              "infer_experiment/$filename"
            else if (filename.indexOf("read_distribution.txt") > 0)             "read_distribution/$filename"
            else if (filename.indexOf("read_duplication.DupRate_plot.pdf") > 0) "read_duplication/$filename"
            else if (filename.indexOf("read_duplication.DupRate_plot.r") > 0)   "read_duplication/rscripts/$filename"
            else if (filename.indexOf("read_duplication.pos.DupRate.xls") > 0)  "read_duplication/dup_pos/$filename"
            else if (filename.indexOf("read_duplication.seq.DupRate.xls") > 0)  "read_duplication/dup_seq/$filename"
            else if (filename.indexOf("inner_distance.txt") > 0)                "inner_distance/$filename"
            else if (filename.indexOf("inner_distance_freq.txt") > 0)           "inner_distance/data/$filename"
            else if (filename.indexOf("inner_distance_plot.r") > 0)             "inner_distance/rscripts/$filename"
            else if (filename.indexOf("inner_distance_plot.pdf") > 0)           "inner_distance/plots/$filename"
            else if (filename.indexOf("junction_plot.r") > 0)                   "junction_annotation/rscripts/$filename"
            else if (filename.indexOf("junction.xls") > 0)                      "junction_annotation/data/$filename"
            else if (filename.indexOf("splice_events.pdf") > 0)                 "junction_annotation/events/$filename"
            else if (filename.indexOf("splice_junction.pdf") > 0)               "junction_annotation/junctions/$filename"
            else if (filename.indexOf("junctionSaturation_plot.pdf") > 0)       "junction_saturation/$filename"
            else if (filename.indexOf("junctionSaturation_plot.r") > 0)         "junction_saturation/rscripts/$filename"
            else filename
        }    
    when:
    !params.skip_qc && !params.skip_rseqc

    input:
    file bam_file from rseqc_bam
    file bed12 from bed12file.collect()

    output:
    file "*" into rseqc_results
    file "*.bam_stat.txt" into bam_stat_for_normlization

    script:
    output = bam_file[0].toString() - ~/(\.bam)?$/
    """    
    infer_experiment.py -i ${bam_file[0]} -r ${bed12} > ${output}.infer_experiment.txt
    junction_annotation.py -i ${bam_file[0]}  -o ${output}.rseqc -r ${bed12}
    bam_stat.py -i ${bam_file[0]}  > ${output}.bam_stat.txt
    junction_saturation.py -i ${bam_file[0]}  -o ${output}.rseqc -r ${bed12} 2> ${output}.junction_annotation_log.txt
    inner_distance.py -i ${bam_file[0]}  -o ${output}.rseqc -r ${bed12}
    read_distribution.py -i ${bam_file[0]}  -r ${bed12} > ${output}.read_distribution.txt
    read_duplication.py -i ${bam_file[0]}  -o ${output}.read_duplication
    """
}


process CreateBedgraph{
    tag "$output"
    publishDir "${params.outdir}/QC/rseqc/", mode: 'link', overwrite: true ,
        saveAs: {filename ->
            if (filename.indexOf("bedgraph") > 0) "bedgraph/$filename"
        }

    input:
    file bam_file from bedgraph_bam

    output:
    file "*.bedgraph" into bedgraph_for_genebody, bedgraph_for_igv

    when:
    !params.skip_createbedgraph

    script:
    output = bam_file[0].toString() - ~/(\.bam)?$/
    """
    bamCoverage -b ${bam_file[0]} --outFileFormat bedgraph -o ${output}.bedgraph -p ${task.cpus}
    """
}

Channel
    .from()
    .concat( tophat2_log, hisat2_log, star_log, fastp_results, fastqc_results, rseqc_results)
    .into{ arranged_qc; qc_results_for_report }

process multiqc{
    publishDir "${params.outdir}/Report/QCReadsReport" , mode: 'link', overwrite: true
    
    when:
    !params.skip_qc

    input:
    file arranged_qc from arranged_qc.collect()

    output:
    file "multiqc*" into multiqc_results

    script:
    """
    multiqc -n multiqc_report_$aligner .
    """
}

index_peakCallingbygroup = params.peakCalling_mode == "group" ? 0 : 1
sorted_bam.groupTuple(by: [0,1]).groupTuple(by: index_peakCallingbygroup)
.map{group,sample,bam,bai -> [group,sample, bam.flatten().sort{
    o1, o2 -> 
    sub_o1 = o1.toString();sub_o1 = sub_o1.substring(sub_o1.lastIndexOf('/') + 1);
    sub_o2 = o2.toString();sub_o2 = sub_o2.substring(sub_o2.lastIndexOf('/') + 1);
    return sub_o1.compareTo(sub_o2);
}, bai.flatten()]}
.into{ macs2_bam; metpeak_bam; matk_bam; meyer_bam}


/*
========================================================================================
                        Step 4 Differential expression analysis
========================================================================================
*/
process FeatureCount{
    label 'analysis'
    publishDir "${params.outdir}/expressionAnalysis/featurecounts", mode: 'link', overwrite: true

    input:
    file bam_bai_file from feacount_bam.collect()
    file formatted_designfile from formatted_designfile.collect()
    file gtf

    output:
    //set file("*input*.count"), val(sample_name) into count_input
    file "*input*.count" into htseq_count_input, htseq_count_input_to_arrange
    file "expression*.matrix" into htseq_results
    file "expression*.count.matrix" into htseq_diffm6a_results

    script:
    println LikeletUtils.print_purple("Generate gene expression matrix by featureCounts and Rscript")
    strand_info = params.stranded == "no" ? "0" : params.stranded == "reverse" ? "2" : "1"
    reads_type = params.single_end ? "" : "-p "
    minmqs = params.featurecount_minMQS
    align_tool = params.aligners
    // """
    // bash $baseDir/bin/htseq_count.sh $gtf $strand_info ${task.cpus}
    // Rscript $baseDir/bin/get_htseq_matrix.R $formatted_designfile ${task.cpus} 
    // """
    """
    for bam_file in *.input*.bam
    do
        featureCounts $reads_type-Q $minmqs -T ${task.cpus} -s $strand_info -a ${gtf} -o \${bam_file/%_sort*/}.txt \${bam_file};
    done
    Rscript $baseDir/bin/generate_featurecount_mat.R $formatted_designfile ${task.cpus} $gtf $align_tool
    """ 
}

process DESeq2{
    label 'analysis'
    tag "$compare_str"

    publishDir "${params.outdir}/expressionAnalysis/DESeq2", mode: 'link', overwrite: true

    input:
    file reads_count_input from htseq_count_input.collect()
    file formatted_designfile from formatted_designfile.collect()
    val compare_str from compareLines_for_DESeq2

    output:
    file "DESeq2*.csv" into deseq2_results
    
    when:
    !params.skip_deseq2 && !params.skip_expression && params.comparefile
    
    script:
    println LikeletUtils.print_purple("Differential expression analysis performed by DESeq2 ($compare_str)")
    """
    Rscript $baseDir/bin/DESeq2.R $formatted_designfile $compare_str $params.aligners
    """ 
}

process EdgeR{
    label 'analysis'
    tag "$compare_str"
    publishDir "${params.outdir}/expressionAnalysis/edgeR", mode: 'link', overwrite: true

    input:
    file reads_count_input from htseq_count_input.collect()
    file formatted_designfile from formatted_designfile.collect()
    val compare_str from compareLines_for_edgeR

    output:
    file "edgeR*.csv" into edgeR_results
    
    when:
    !params.skip_edger && !params.skip_expression && params.comparefile

    script:
    println LikeletUtils.print_purple("Differential expression analysis performed by EdgeR ($compare_str)")
    """
    Rscript $baseDir/bin/edgeR.R $formatted_designfile $compare_str $params.aligners
    """ 
}


/*
========================================================================================
                            Step 5 Peak Calling
========================================================================================
*/ 
/*
 * STEP 5 - 1  Peak Calling------MetPeak, MACS2, MATK
*/
process Metpeak {
    tag "$peakcalling_tag"
    label 'onecore_peak'
    publishDir "${params.outdir}/peakCalling/metpeak", mode: 'link', overwrite: true

    input:
    set val(group), val(sample_id), file(bam_vector), file(bai_vector) from metpeak_bam
    file gtf

    output:
    file "*" into metpeak_results
    file "metpeak*_normalized.bed" into metpeak_nomarlized_bed

    when:
    !params.skip_metpeak && !params.skip_peakCalling

    script:
    input = []; ip = []
    bam_vector.eachWithIndex{val, ix -> (ix & 1 ? ip : input) << val}
    input_bam = input.join(','); ip_bam = ip.join(',')
    peakcalling_tag = params.peakCalling_mode == "group" ? "group_" + group : sample_id
    """
    echo $input_bam
    Rscript $baseDir/bin/MeTPeak.R $input_bam $ip_bam $peakcalling_tag $gtf
    """
}

process Macs2{
    tag "$peakcalling_tag"
    label 'onecore_peak'
    publishDir "${params.outdir}/peakCalling/macs2", mode: 'link', overwrite: true

    input:
    file fasta
    set val(group), val(sample_id), file(bam_vector), file(bai_vector) from macs2_bam

    output:
    file "macs2*.{xls,narrowPeak,summits}" into macs2_results
    file "macs2*_normalized.bed" into macs2_nomarlized_bed

    when:
    !params.skip_macs2 && !params.skip_peakCalling

    script:
    input = []; ip = []
    bam_vector.eachWithIndex{val, ix -> (ix & 1 ? ip : input) << val}
    input_bam = input.join(' '); ip_bam = ip.join(' ')
    input_file_count = input.size()
    ip_file_count = ip.size()
    peakcalling_tag = params.peakCalling_mode == "group" ? "group_" + group : sample_id
    arguments = params.peak_threshold == "high" ? "-p 1e-6 --keep-dup 5" : params.peak_threshold == "medium" ? "-q 0.01 --keep-dup 5" : "-q 0.05 --keep-dup 3"
    """
    genome_size=\$(faCount $fasta | tail -1 | awk '{print \$2-\$7}')
    if [ $ip_file_count -gt 1 ]; then samtools merge -f macs2_${peakcalling_tag}_ip.bam $ip_bam; fi
    if [ $input_file_count -gt 1 ]; then samtools merge -f macs2_${peakcalling_tag}_input.bam $input_bam; fi
    macs2 callpeak -t *${peakcalling_tag}*ip*.bam -c *${peakcalling_tag}*input*.bam -g \$genome_size -n macs2_${peakcalling_tag} $arguments -f BAM --nomodel
    awk -v OFS="\\t" '{print \$1,\$2,\$3,\$1":"\$2"-"\$3,10^-\$8}' macs2_${peakcalling_tag}_peaks.narrowPeak > macs2_${peakcalling_tag}_normalized.bed
    mv macs2_${peakcalling_tag}_summits.bed macs2_${peakcalling_tag}.summits
    """ 
}


process MATKpeakCalling {
    tag "$peakcalling_tag"
    label 'onecore_peak'
    publishDir "${params.outdir}/peakCalling/MATK", mode: 'link', overwrite: true

    input:
    set val(group), val(sample_id), file(bam_vector), file(bai_vector) from matk_bam
    file gtf

    output:
    file "*" into matk_results
    file "MATK*_normalized.bed" into matk_nomarlized_bed

    when:
    !params.skip_matk && !params.skip_peakCalling

    script:
    matk_jar = params.matk_jar
    input = []; ip = []
    bam_vector.eachWithIndex{val, ix -> (ix & 1 ? ip : input) << val}
    input_bam = input.join(';'); ip_bam = ip.join(';')
    input_file_count = input.size()
    peakcalling_tag = params.peakCalling_mode == "group" ? "group_" + group : sample_id
    arguments = params.peak_threshold == "high" ? "-q 0.01" : params.peak_threshold == "medium" ? "-q 0.05" : "-q 0.1"
    """
    #export OMP_NUM_THREADS=${task.cpus}
    if [ ! -f "$matk_jar" ]; then
        echo "Cannot find matk.jar. Please check the param of matk_jar" 1>&2
        exit 1
    fi
    java -jar $matk_jar -peakCalling $arguments -c $input_file_count -ip "$ip_bam" -input "$input_bam" -out MATK_${peakcalling_tag}.bed
    awk 'BEGIN{FS="\\t";OFS="\\t"}{print \$1,\$2,\$3,\$1":"\$2"-"\$3,\$5}' MATK_${peakcalling_tag}.bed > MATK_${peakcalling_tag}_normalized.bed
    """    
}

/*
 * PREPROCESSING - Create chromsizesfile for meyer
 * NEED genome.fa
 */
if( params.fasta ){
    process MeyerPrepration{
        label 'build_index'
        tag "onecore_peak"
        publishDir path: { params.saveReference ? "${params.outdir}/Genome/meyerPrepration" : params.outdir },
                saveAs: { params.saveReference ? it : null }, mode: 'copy'       

        input:
        file fasta
        file chromsizesfile from chromsizesfile.collect()

        output:
        file "chrName.txt" into chrNamefile
        file "genome.bin25.bed" into bin25file
        file "genomebin" into genomebin

        when:
        !params.skip_meyer && !params.skip_peakCalling

        shell:
        '''
        awk '{print $1}' !{chromsizesfile} > chrName.txt
        mkdir genomebin
        bedtools makewindows -g !{chromsizesfile} -w 25 > genome.bin25.bed
        awk '{print "cat genome.bin25.bed | grep "$1" > genomebin/"$1".bin25.bed"}' chrName.txt | xargs -iCMD -P!{task.cpus} bash -c CMD
        '''
    }
}else {
    exit 1, println LikeletUtils.print_red("Chromsizes file cannot build due to lack of " + params.fasta)
}
process Meyer{
    tag "$peakcalling_tag"
    label 'peak_calling'
    publishDir "${params.outdir}/peakCalling/meyer", mode: 'link', overwrite: true

    input:
    set val(group), val(sample_id), file(bam_vector), file(bai_vector) from meyer_bam
    file chrNamefile from chrNamefile
    file bin25file from bin25file
    file genomebin from genomebin

    output:
    file "meyer*.bed" into meyer_results
    file "meyer*_normalized.bed" into meyer_nomarlized_bed

    when:
    !params.skip_meyer && !params.skip_peakCalling

    shell:
    
    // if( flag_peakCallingbygroup ){
    //     println LikeletUtils.print_purple("Peak Calling performed by Meyer in group mode")
    // }else{
    //     println LikeletUtils.print_purple("Peak Calling performed by Meyer in independent mode")
    // }
    input = []; ip = []
    bam_vector.eachWithIndex{val, ix -> (ix & 1 ? ip : input) << val}
    input_bam = input.join(' '); ip_bam = ip.join(' ')
    input_file_count = input.size()
    ip_file_count = ip.size()
    peakcalling_tag = params.peakCalling_mode == "group" ? "group_" + group : sample_id
    arguments = params.peak_threshold == "high" ? "-q 0.01" : params.peak_threshold == "medium" ? "-q 0.05" : "-q 0.1"
    '''
    cp !{baseDir}/bin/meyer.py ./
    if [ !{ip_file_count} -gt 1 ]; then samtools merge -f meyer_!{peakcalling_tag}_ip.bam !{ip_bam}; fi
    if [ !{input_file_count} -gt 1 ]; then samtools merge -f meyer_!{peakcalling_tag}_input.bam !{input_bam}; fi    
    input_bam=*!{peakcalling_tag}*input*.bam
    ip_bam=*!{peakcalling_tag}*ip*.bam
    prefix=meyer_!{peakcalling_tag}
    genomebin_dir="genomebin/"
    peak_windows_number=$(wc -l !{bin25file}| cut -d " " -f 1)
    input_total_reads_count=$(samtools view -c $input_bam)
    ip_total_reads_count=$(samtools view -c $ip_bam)
    mkdir $prefix.tmp "$prefix.tmp/input" "$prefix.tmp/ip"
    awk -v bam="$input_bam" -v pre="$prefix" '
        {print "samtools view -b "bam" "$1 ">./"pre".tmp/input/"$1".bam; \
        bamToBed -split -i < ./"pre".tmp/input/"$1".bam>./"pre".tmp/input/"$1".bed; \
        sortBed -i ./"pre".tmp/input/"$1".bed | intersectBed  -a '${genomebin_dir}'"$1".bin25.bed -b - -sorted -c > ./"pre".tmp/input/"$1".bin25.txt"}' !{chrNamefile} \
        | xargs -iCMD -P!{task.cpus} bash -c CMD
    awk -v bam="$ip_bam" -v pre="$prefix"  '
        {print "samtools view -b "bam" "$1 ">./"pre".tmp/ip/"$1".bam; \
        bamToBed -split -i < ./"pre".tmp/ip/"$1".bam>./"pre".tmp/ip/"$1".bed; \
        sortBed -i ./"pre".tmp/ip/"$1".bed | intersectBed  -a '${genomebin_dir}'"$1".bin25.bed -b - -sorted -c > ./"pre".tmp/ip/"$1".bin25.txt"}' !{chrNamefile} \
        | xargs -iCMD -P!{task.cpus} bash -c CMD
    echo "cal pval for each 25bp bin"
    awk -v pre="$prefix" '
    {print "python meyer.py ./"pre".tmp/input/"$1".bin25.txt ./"pre".tmp/ip/"$1".bin25.txt '$input_total_reads_count' '$ip_total_reads_count' '$peak_windows_number' ./"pre".tmp/ip/"$1".m6A.meyer.pval.txt"}' !{chrNamefile} \
    |xargs -iCMD -P!{task.cpus} bash -c CMD
    cat $prefix.tmp/ip/*.m6A.meyer.pval.txt > ${prefix}.bed
    awk 'BEGIN{FS="\t";OFS="\t"}{print $1,$2,$3,$1":"$2"-"$3,$4}' ${prefix}.bed > ${prefix}_normalized.bed
    rm -rf $prefix.tmp
    ''' 
}

/*
========================================================================================
                        Step 6 Merge Peak AND Analysis
========================================================================================
*/
/*
 * STEP 6-1 Merge Peak
*/

Channel
    .from()
    .concat(metpeak_nomarlized_bed, macs2_nomarlized_bed, matk_nomarlized_bed, meyer_nomarlized_bed)
    .into {merged_bed ; bed_for_annotation}

process PeakMerge {
    label 'analysis'
    publishDir "${params.outdir}/peakCalling/mergedBed", mode: 'link', overwrite: true,
        saveAs: {filename ->
            if (filename.indexOf("bed") > 0) "$filename"
        }
    
    input:
    file peak_bed from merged_bed.collect()
    file formatted_designfile from formatted_designfile.collect()

    output:
    file "*merged*.bed" into merge_result
    file "*merged_group*.bed" into group_merged_bed
    file "*_merged_allpeaks.bed" into all_merged_bed

    script:
    flag_peakCallingbygroup = params.peakCalling_mode == "group" ? 1 : 0
    peakCalling_tools_count = (params.skip_metpeak ? 0 : 1).toInteger() + (params.skip_macs2 ? 0 : 1).toInteger() + (params.skip_matk ? 0 : 1).toInteger() + (params.skip_meyer ? 0 : 1).toInteger()
    peakMerged_mode = params.peakMerged_mode
    if ( peakMerged_mode == "rank" )  
        println LikeletUtils.print_purple("Start merge peaks by RobustRankAggreg")
    else if ( peakMerged_mode == "mspc" )  
        println LikeletUtils.print_purple("Start merge peaks by MSPC")
    else
        println LikeletUtils.print_purple("Start merge peaks by " + peakMerged_mode )

    """
    cp ${baseDir}/bin/normalize_peaks.py ./
    if [ ${peakMerged_mode} == "rank" ]; then 
        cp $baseDir/bin/merge_peaks_by_rank.R ./
        bash $baseDir/bin/merge_peaks_by_rank.sh $formatted_designfile ${task.cpus} $flag_peakCallingbygroup $peakCalling_tools_count
    elif [ ${peakMerged_mode} == "mspc" ]; then
        bash $baseDir/bin/merge_peaks_by_mspc.sh $formatted_designfile ${task.cpus} $flag_peakCallingbygroup $peakCalling_tools_count mspc_results
    elif [ ${peakMerged_mode} == "macs2" ]||[ ${peakMerged_mode} == "metpeak" ]||[ ${peakMerged_mode} == "MATK" ]||[ ${peakMerged_mode} == "meyer" ]; then 
        bash $baseDir/bin/merge_peaks_by_bedtools.sh $formatted_designfile ${task.cpus} $flag_peakCallingbygroup $peakCalling_tools_count $peakMerged_mode
    else
        echo -e "Please check your value of peakMerged_mode: $peakMerged_mode"
    fi
    whether_nopeaks=\$(wc -l *merged*.bed | awk '\$1==0{print "error"}' | uniq)
    if [[ \$whether_nopeaks == "error" ]] ;then 
        echo "There is no peaks in one of the merged peaks files" 1>&2
        echo "Merge Peaks by "${peakMerged_mode}" may not be suitable for your data." 1>&2
        exit 1
    fi
    """
}

Channel
    .from()
    .mix( group_merged_bed, all_merged_bed )
    .flatten()
    .into{ annotate_collection; motif_collection; bed_collect_for_arrange_results}
annotate_collection.mix(bed_for_annotation).toList().flatten().set{beds_anno}

process BedAnnotated{
    tag "${all_bed.baseName}"
    label 'onecore_peak'
    publishDir "${params.outdir}/m6AAnalysis/AnnotatedPeaks", mode: 'link', overwrite: true
    
    input:
    file all_bed from beds_anno
    file fasta
    file gtf

    output:
    file "annotatedbygtf/*" into annotation_results_xy, anno_for_quantification, anno_for_diffreport
    file "*_annotatedbyhomer.bed" into annotation_results_homer
    
    when:
    !params.skip_annotation

    script:
    annotated_script_dir = baseDir + "/bin"
    bed_prefix = all_bed.toString() - ~/(\.bed)?$/
    //Skip Peak Calling Tools Setting
    """
    # Annotation Peaks
    cp ${annotated_script_dir}/intersec.pl ./
    cp ${annotated_script_dir}/m6A_annotate_forGTF_xingyang2.pl ./
    mkdir annotatedbygtf
    perl m6A_annotate_forGTF_xingyang2.pl ${gtf} ${all_bed} annotatedbygtf/${bed_prefix} 
    annotatePeaks.pl ${all_bed} ${fasta} -gtf ${gtf} > ${bed_prefix}_annotatedbyhomer.bed
    """
}

process MotifSearching {
    label 'onecore_peak'
    tag "${bed_file.baseName}"
    publishDir "${params.outdir}/m6AAnalysis/motif", mode: 'link', overwrite: true
    
    input:
    file bed_file from motif_collection
    file chromsizesfile from chromsizesfile.collect()
    file bed12 from bed12file.collect()
    file fasta
    file gtf

    output:
    file "*_{dreme,homer}" into motif_results, motif_results_for_report

    when:
    !params.skip_motif

    script:
    motif_file_dir = baseDir + "/bin"
    bed_prefix = bed_file.baseName
    length = params.motiflength
    println LikeletUtils.print_purple("Motif analysis is going on by Homer")
    """
    # cp ${motif_file_dir}/m6A_motif.meme ./
    sort -k5,5 -g ${bed_file} | awk 'FNR <= 2000{ print \$1"\\t"\$2"\\t"\$3}' > ${bed_prefix}.location
    intersectBed -wo -a ${bed_prefix}.location -b $gtf | awk -v OFS="\\t" '{print \$1,\$2,\$3,"*","*",\$10}' | sort -k1,2 | uniq > ${bed_prefix}_bestpeaks.bed
    fastaFromBed -name+ -split -s -fi $fasta -bed ${bed_prefix}_bestpeaks.bed > ${bed_prefix}_bestpeaks.fa
    # ame -oc ${bed_prefix}_ame ${bed_prefix}_bestpeaks.fa m6A_motif.meme
    shuffleBed -incl ${bed12} -seed 12345 -noOverlapping -i ${bed_prefix}_bestpeaks.bed -g ${chromsizesfile} > ${bed_prefix}_random_peak.bed
    fastaFromBed -name+ -split -s -fi $fasta -bed ${bed_prefix}_random_peak.bed > ${bed_prefix}_random_peak.fa
    findMotifs.pl ${bed_prefix}_bestpeaks.fa fasta ${bed_prefix}_homer -fasta ${bed_prefix}_random_peak.fa -p ${task.cpus} \
        -len $length -S 10 -rna -dumpFasta > ${bed_prefix}_homer_run.log 2>&1
    """
}

process PeaksMotifReport {
    publishDir "${params.outdir}/Report/PeaksMotifReport", mode: 'link', overwrite: true
    
    input:
    file motif from motif_results.collect()
    file annotation_files from annotation_results_xy.collect()
    file formatted_designfile from formatted_designfile.collect()

    output:
    file "*.{html,pdf}" into peaksMotifReport

    script:
    peakMerged_mode = params.peakMerged_mode 
    peakCalling_mode = params.peakCalling_mode
    qcPeaksRData = "PeakMotifPlot.RData"
    """
    cp $baseDir/bin/Peaks_Motif_Report.rmd ./
    Rscript $baseDir/bin/Peaks_Motif_Report.R $formatted_designfile $peakMerged_mode $peakCalling_mode $qcPeaksRData
    R -e "load(\\"$qcPeaksRData\\");rmarkdown::render('Peaks_Motif_Report.rmd',output_file='Peaks_Motif_Report_${peakMerged_mode}.html')"
    """
}

anno_for_quantification.flatten().filter( ~/.*merged_allpeaks.anno.txt$/).set{methylation_annotaion_file}
process PeaksQuantification{
    label 'analysis'
    publishDir "${params.outdir}/m6AAnalysis/m6AQuantification", mode: 'link', overwrite: true
    
    input:
    file merged_bed from all_merged_bed.collect()
    file htseq_count_file from htseq_count_input.collect()
    file bam_bai_file from peakquan_bam.collect()
    file formatted_designfile from formatted_designfile.collect()
    file annotation_file from methylation_annotaion_file.collect()
    file gtf

    output:
    file "*.{matrix,count}" into quantification_results, quantification_matrix

    when:
    !params.skip_peakCalling

    script:
    matk_jar = params.matk_jar
    methylation_analysis_mode = params.methylation_analysis_mode
    if ( methylation_analysis_mode == "Wilcox-test" )  
        println LikeletUtils.print_purple("Generate m6A quantification matrix by bedtools")
    else if ( methylation_analysis_mode == "QNB" )  
        println LikeletUtils.print_purple("Generate m6A quantification matrix by QNB")
    else if ( methylation_analysis_mode == "MATK" )
        println LikeletUtils.print_purple("Generate m6A quantification matrix by MATK")
    else if ( methylation_analysis_mode == "edgeR" )  
        println LikeletUtils.print_purple("Generate m6A quantification matrix by bedtools")
    else if ( methylation_analysis_mode == "DESeq2" )
        println LikeletUtils.print_purple("Generate m6A quantification matrix by bedtools")
    """
    if [ ${methylation_analysis_mode} == "DESeq2" ]||[ ${methylation_analysis_mode} == "edgeR" ]; then 
        # PeaksQuantification by LRT
        bash $baseDir/bin/bed_count.sh ${formatted_designfile} ${task.cpus} ${merged_bed} bam_stat_summary.txt
        Rscript $baseDir/bin/bedtools_quantification.R $formatted_designfile bam_stat_summary.txt
    else
        case ${methylation_analysis_mode} in 
        Wilcox-test)
            bash $baseDir/bin/bed_count.sh ${formatted_designfile} ${task.cpus} ${merged_bed} bam_stat_summary.txt
            Rscript $baseDir/bin/bedtools_quantification.R $formatted_designfile bam_stat_summary.txt
            ;;
        QNB)
            bash $baseDir/bin/bed_count.sh ${formatted_designfile} ${task.cpus} ${merged_bed} bam_stat_summary.txt
            Rscript $baseDir/bin/QNB_quantification.R $formatted_designfile ${task.cpus}
            ;;
        MATK)
            export OMP_NUM_THREADS=${task.cpus}
            if [ ! -f "$matk_jar" ]; then
                echo "Cannot find matk.jar. Please check the param of matk_jar" 1>&2
                exit 1
            fi
            bash $baseDir/bin/MATK_quantification.sh $matk_jar $gtf $formatted_designfile ${merged_bed} 1
            ;;
        *)
            echo ${methylation_analysis_mode}" is not Wilcox-test, QNB, MATK, DESeq2 or edgeR"
            exit 1
            ;;
        esac
    fi
    head -1 *_quantification.matrix |sed 's/^\\t//'  |awk -F "\\t" '{print "ID\\tGene_symbol\\t"\$0}' > tmp.header.file
    sed '1d' *_quantification.matrix | sort > tmp.quantification.file
    awk 'BEGIN{FS="\\t";OFS="\\t"}{print \$4,\$15,\$11}' ${annotation_file} | sort | join -t \$'\t' -e 'NA' -a1 -o 1.1 -o 2.2 -o 2.3 tmp.quantification.file - >  tmp.annotation.file
    join -t \$'\t' tmp.annotation.file tmp.quantification.file | cat tmp.header.file - > *_quantification.matrix 
    """
}


process diffm6APeak{
    label 'analysis'
    tag "$compare_str"
    publishDir "${params.outdir}/m6AAnalysis/diffm6A", mode: 'link', overwrite: true
    
    input:
    //file peak_bed from group_merged_bed.collect()
    file merged_bed from all_merged_bed.collect()
    file bam_bai_file from diffpeak_bam.collect()
    file formatted_designfile from formatted_designfile.collect()
    file count_matrix from quantification_matrix.collect()
    file exp_matrix from htseq_diffm6a_results.collect()
    file gtf
    val compare_str from compareLines_for_diffm6A

    output:
    file "*diffm6A*.txt" into diffm6A_results

    when:
    !params.skip_diffpeakCalling && params.comparefile

    script:
    matk_jar = params.matk_jar
    methylation_analysis_mode = params.methylation_analysis_mode
    if ( methylation_analysis_mode == "Wilcox-test" )  
        println LikeletUtils.print_purple("Differential m6A analysis is going on by bedtools")
    else if ( methylation_analysis_mode == "QNB" )  
        println LikeletUtils.print_purple("Differential m6A analysis is going on by QNB")
    else if ( methylation_analysis_mode == "MATK" )
        println LikeletUtils.print_purple("Differential m6A analysis is going on by MATK")
    else if ( methylation_analysis_mode == "edgeR" )  
        println LikeletUtils.print_purple("Differential m6A analysis is going on by edgeR")
    else if ( methylation_analysis_mode == "DESeq2" )
        println LikeletUtils.print_purple("Differential m6A analysis is going on by DESeq2")
    """
    case ${methylation_analysis_mode} in 
        Wilcox-test)
            Rscript $baseDir/bin/bedtools_diffm6A.R $formatted_designfile bedtools_quantification.matrix $compare_str
            ;;
        QNB)
            Rscript $baseDir/bin/QNB_diffm6A.R $formatted_designfile ${merged_bed} $compare_str   
            ;;
        MATK)
            export OMP_NUM_THREADS=${task.cpus}
            if [ ! -f "$matk_jar" ]; then
                echo "Cannot find matk.jar. Please check the param of matk_jar" 1>&2
                exit 1
            fi
            bash $baseDir/bin/MATK_diffm6A.sh $matk_jar $formatted_designfile $gtf $compare_str $merged_bed
            ;;
        edgeR)
            Rscript $baseDir/bin/GLM_edgeR_DM.R $formatted_designfile $compare_str bedtools_quantification.matrix $exp_matrix   
            ;;
        DESeq2)
            Rscript $baseDir/bin/GLM_DESeq2_DM.R $formatted_designfile $compare_str ${task.cpus} bedtools_quantification.matrix $exp_matrix 
            ;;
        *)
            echo ${methylation_analysis_mode}" is not Wilcox-test, QNB, MATK, DESeq2 or edgeR"
            exit 1
            ;;
    esac   
    """ 
}

process SingleNucleotidePrediction{
    label 'analysis'
    publishDir "${params.outdir}/m6AAnalysis/m6APredictionSites", mode: 'link', overwrite: true
    
    input:
    file peak_bed from group_merged_bed.collect()
    file group_bed from all_merged_bed.collect()
    file formatted_designfile from formatted_designfile.collect()
    file bam_bai_file from sng_bam.collect()
    file fasta
    file gtf

    output:
    file "m6A_sites*.bed" into prediction_results

    when:
    !params.skip_m6Aprediction

    script:
    matk_jar = params.matk_jar
    println LikeletUtils.print_purple("SignleNucleotide Prediction analysis is going on by MATK")
    """
    export OMP_NUM_THREADS=${task.cpus}
    bash $baseDir/bin/m6Aprediction.sh $matk_jar $formatted_designfile $fasta $gtf
    """
}


Channel
    .from()
    .concat( quantification_results, motif_results_for_report, diffm6A_results, 
        htseq_count_input_to_arrange, 
        anno_for_diffreport, prediction_results, bed_collect_for_arrange_results,
        multiqc_results, deseq2_results, edgeR_results
    )
    .set{ results_arrange }

process DiffReport {
    publishDir "${params.outdir}/Report" , mode: 'link', overwrite: true,
        saveAs: {filename ->
                 if (filename.indexOf(".html") > 0)  "diffReport/$filename"
                 else if (filename.indexOf(".pdf") > 0)  "diffReport/$filename"
                 else "ReportRData/$filename"
        }        
    input:
    file results from results_arrange.collect()
    file formatted_designfile from formatted_designfile.collect()
    val compare_info from compareLines_for_arranged_result.collect()
    
    output:
    file "*.m6APipe" into m6APipe_result
    file "*.{html,pdf}" into diffReport_result

    when:
    !params.skip_annotation && !params.skip_expression && !params.skip_diffpeakCalling && !params.skip_peakCalling && params.comparefile
    
    script:
    methylation_analysis_mode = params.methylation_analysis_mode
    expression_analysis_mode = params.expression_analysis_mode
    peakMerged_mode = params.peakMerged_mode
    diffReportRData = "DiffReport.RData"
    """
    cp $baseDir/bin/DiffReport.rmd ./
    if [ "$compare_info" != "[two_groups]" ]; then
        echo $compare_info | sed 's/^\\[//g' | sed 's/\\]\$//g' | sed s/[[:space:]]//g > compare_info
    else
        echo \$(awk 'BEGIN{FS=","}NR>1{print \$4}' $formatted_designfile |sort|uniq|awk 'NR==1{printf \$0"_vs_"}NR==2{print \$0}') > compare_info
    fi
    Rscript $baseDir/bin/arranged_results.R $formatted_designfile compare_info $methylation_analysis_mode $expression_analysis_mode $peakMerged_mode
    Rscript $baseDir/bin/DiffReport.R *.m6APipe $diffReportRData
    R -e "load(\\"$diffReportRData\\");rmarkdown::render('DiffReport.rmd',output_file='DiffReport_${peakMerged_mode}_${methylation_analysis_mode}_${expression_analysis_mode}.html')"
    rm Rplots.pdf
    """
}

process CreateIGVjs {
    publishDir "${params.outdir}/Report" , mode: 'link', overwrite: true,
        saveAs: {filename ->
                 if (filename.indexOf(".html") > 0)  "Igv_js/$filename"
                 else if (filename.indexOf(".pdf") > 0)  "Igv_js/$filename"
                 else "Igv_js/$filename"
        }        
    input:
    file m6APipe_result from m6APipe_result
    file fasta 
    file gtf
    file formatted_designfile from formatted_designfile.collect()
    file group_bed from group_merged_bed.collect()
    file all_bed from all_merged_bed.collect()
    file bedgraph from bedgraph_for_igv.collect()
    
    output:
    file "*" into igv_js

    script:    
    igv_fasta = fasta.baseName.toString() + ".igv.fa"
    igv_gtf = gtf.baseName.toString() + ".igv.gtf"
    merged_allpeaks_igvfile = all_bed.baseName.toString() + ".igv.bed"
    """
    ls -l $fasta | awk -F "> " '{print "ln -s "\$2" ./'$igv_fasta'"}' | bash
    ls -l $gtf | awk -F "> " '{print "ln -s "\$2" ./'$igv_gtf'"}' | bash
    ls -l $m6APipe_result | awk '{print "ln -s "\$11" initial.m6APipe"}' | bash
    ls -l $group_bed $all_bed | awk '{sub(".bed\$",".igv.bed",\$9);print "ln -s "\$11,\$9}' | bash
    ls -l $bedgraph | awk '{sub(".bedgraph\$",".igv.bedgraph",\$9);print "ln -s "\$11,\$9}' | bash
    samtools faidx $igv_fasta
    bash $baseDir/bin/create_IGV_js.sh $igv_fasta $igv_gtf $merged_allpeaks_igvfile $formatted_designfile
    """
}

/*
 * Parse software version numbers
 */
process get_software_versions {
    publishDir "${params.outdir}/pipeline_info", mode: 'copy',
        saveAs: { filename ->
            if (filename.indexOf(".csv") > 0) filename
            else null
        }

    output:
    file 'software_versions_mqc.yaml' into software_versions_yaml
    file "software_versions.csv"

    script:
    // TODO nf-core: Get all tools to print their version number here
    """
    echo $workflow.manifest.version > v_pipeline.txt
    echo $workflow.nextflow.version > v_nextflow.txt
    fastqc --version > v_fastqc.txt
    multiqc --version > v_multiqc.txt
    python ${baseDir}/bin/scrape_software_versions.py &> software_versions_mqc.yaml
    """
}



/*
 * Completion e-mail notification
 */
def summary = [:]
if (workflow.revision) summary['Pipeline Release'] = workflow.revision
summary['Run Name']         = custom_runName ?: workflow.runName
// TODO nf-core: Report custom parameters here
summary['Reads']            = params.input
summary['Fasta Ref']        = params.fasta
// summary['Data Type']        = params.single_end ? 'Single-End' : 'Paired-End'
summary['Max Resources']    = "$params.max_memory memory, $params.max_cpus cpus, $params.max_time time per job"
if (workflow.containerEngine) summary['Container'] = "$workflow.containerEngine - $workflow.container"
summary['Output dir']       = params.outdir
summary['Launch dir']       = workflow.launchDir
summary['Working dir']      = workflow.workDir
summary['Script dir']       = workflow.projectDir
summary['User']             = workflow.userName
if (params.config_profile_description) summary['Config Description'] = params.config_profile_description
if (params.config_profile_contact)     summary['Config Contact']     = params.config_profile_contact
if (params.config_profile_url)         summary['Config URL']         = params.config_profile_url
if (params.email || params.email_on_fail) {
  summary['E-mail Address']    = params.email
  summary['E-mail on failure'] = params.email_on_fail
  summary['MultiQC maxsize']   = params.maxMultiqcEmailFileSize
}
workflow.onComplete {

    // Set up the e-mail variables
    def subject = "[MeRIPseqPpipe] Successful: $workflow.runName"
    if (!workflow.success) {
      subject = "[MeRIPseqPpipe] FAILED: $workflow.runName"
    }
    def email_fields = [:]
    email_fields['version'] = workflow.manifest.version
    email_fields['runName'] = custom_runName ?: workflow.runName
    email_fields['success'] = workflow.success
    email_fields['dateComplete'] = workflow.complete
    email_fields['duration'] = workflow.duration
    email_fields['exitStatus'] = workflow.exitStatus
    email_fields['errorMessage'] = (workflow.errorMessage ?: 'None')
    email_fields['errorReport'] = (workflow.errorReport ?: 'None')
    email_fields['commandLine'] = workflow.commandLine
    email_fields['projectDir'] = workflow.projectDir
    email_fields['summary'] = summary
    email_fields['summary']['Date Started'] = workflow.start
    email_fields['summary']['Date Completed'] = workflow.complete
    email_fields['summary']['Pipeline script file path'] = workflow.scriptFile
    email_fields['summary']['Pipeline script hash ID'] = workflow.scriptId
    if (workflow.repository) email_fields['summary']['Pipeline repository Git URL'] = workflow.repository
    if (workflow.commitId) email_fields['summary']['Pipeline repository Git Commit'] = workflow.commitId
    if (workflow.revision) email_fields['summary']['Pipeline Git branch/tag'] = workflow.revision
    if (workflow.container) email_fields['summary']['Docker image'] = workflow.container
    email_fields['summary']['Nextflow Version'] = workflow.nextflow.version
    email_fields['summary']['Nextflow Build'] = workflow.nextflow.build
    email_fields['summary']['Nextflow Compile Timestamp'] = workflow.nextflow.timestamp

    // TODO nf-core: If not using MultiQC, strip out this code (including params.maxMultiqcEmailFileSize)
    // On success try attach the multiqc report
    def mqc_report = null
    try {
        if (workflow.success) {
            mqc_report = multiqc_results.getVal()
            if (mqc_report.getClass() == ArrayList) {
                log.warn "[MeRIPseqPpipe] Found multiple reports from process 'multiqc', will use only one"
                mqc_report = mqc_report[0]
            }
        }
    } catch (all) {
        log.warn "[MeRIPseqPpipe] Could not attach MultiQC report to summary email"
    }

    // Check if we are only sending emails on failure
    email_address = params.email
    if (!params.email && params.email_on_fail && !workflow.success) {
        email_address = params.email_on_fail
    }

    // Render the TXT template
    def engine = new groovy.text.GStringTemplateEngine()
    def tf = new File("$baseDir/assets/email_template.txt")
    def txt_template = engine.createTemplate(tf).make(email_fields)
    def email_txt = txt_template.toString()

    // Render the HTML template
    def hf = new File("$baseDir/assets/email_template.html")
    def html_template = engine.createTemplate(hf).make(email_fields)
    def email_html = html_template.toString()

    // Render the sendmail template
    def smail_fields = [ email: email_address, subject: subject, email_txt: email_txt, email_html: email_html, baseDir: "$baseDir", mqcFile: mqc_report, mqcMaxSize: params.maxMultiqcEmailFileSize.toBytes() ]
    def sf = new File("$baseDir/assets/sendmail_template.txt")
    def sendmail_template = engine.createTemplate(sf).make(smail_fields)
    def sendmail_html = sendmail_template.toString()

    // Send the HTML e-mail
    if (email_address) {
        try {
          if ( params.plaintext_email ){ throw GroovyException('Send plaintext e-mail, not HTML') }
          // Try to send HTML e-mail using sendmail
          [ 'sendmail', '-t' ].execute() << sendmail_html
          log.info "[MeRIPseqPpipe] Sent summary e-mail to $email_address (sendmail)"
        } catch (all) {
          // Catch failures and try with plaintext
          [ 'mail', '-s', subject, email_address ].execute() << email_txt
          log.info "[MeRIPseqPpipe] Sent summary e-mail to $email_address (mail)"
        }
    }

    // Write summary e-mail HTML to a file
    def output_d = new File( "${params.outdir}/pipeline_info/" )
    if (!output_d.exists()) {
      output_d.mkdirs()
    }
    def output_hf = new File( output_d, "pipeline_report.html" )
    output_hf.withWriter { w -> w << email_html }
    def output_tf = new File( output_d, "pipeline_report.txt" )
    output_tf.withWriter { w -> w << email_txt }

    c_reset = params.monochrome_logs ? '' : "\033[0m";
    c_purple = params.monochrome_logs ? '' : "\033[0;35m";
    c_green = params.monochrome_logs ? '' : "\033[0;32m";
    c_red = params.monochrome_logs ? '' : "\033[0;31m";

    if (workflow.stats.ignoredCount > 0 && workflow.success) {
      log.info "${c_purple}Warning, pipeline completed, but with errored process(es) ${c_reset}"
      log.info "${c_red}Number of ignored errored process(es) : ${workflow.stats.ignoredCount} ${c_reset}"
      log.info "${c_green}Number of successfully ran process(es) : ${workflow.stats.succeedCount} ${c_reset}"
    }

    if (workflow.success) {
        log.info "${c_purple}[MeRIPseqPpipe]${c_green} Pipeline completed successfully${c_reset}"
    } else {
        checkHostname()
        log.info "${c_purple}[MeRIPseqPpipe]${c_red} Pipeline completed with errors${c_reset}"
    }

}

def checkHostname(){
    def c_reset = params.monochrome_logs ? '' : "\033[0m"
    def c_white = params.monochrome_logs ? '' : "\033[0;37m"
    def c_red = params.monochrome_logs ? '' : "\033[1;91m"
    def c_yellow_bold = params.monochrome_logs ? '' : "\033[1;93m"
    if (params.hostnames) {
        def hostname = "hostname".execute().text.trim()
        params.hostnames.each { prof, hnames ->
            hnames.each { hname ->
                if (hostname.contains(hname) && !workflow.profile.contains(prof)) {
                    log.error "====================================================\n" +
                            "  ${c_red}WARNING!${c_reset} You are running with `-profile $workflow.profile`\n" +
                            "  but your machine hostname is ${c_white}'$hostname'${c_reset}\n" +
                            "  ${c_yellow_bold}It's highly recommended that you use `-profile $prof${c_reset}`\n" +
                            "============================================================"
                }
            }
        }
    }
}


