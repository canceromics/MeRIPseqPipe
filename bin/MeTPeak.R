## Rscript MeTPeak.R <designfile> <gtf> <THREAD_NUM> <flag_peakCallingbygroup> eg. Rscript MeTPeak.R designfile.txt genes.gtf 10
### designfile: Sample_id, Input_filename, IP_filename, group.id
### flag_peakCallingbygroup: 1(group) 0(sample)
library(MeTPeak)
library(parallel)
args <- commandArgs(T)
input.bam.vec <- unlist(strsplit(args[1], split=','))
ip.bam.vec <- unlist(strsplit(args[2], split=','))
group.id <- args[3]
gtf <- args[4]

##Running MeTPeak and rename the output name
metpeak(GENE_ANNO_GTF = gtf,
        IP_BAM = ip.bam.vec,
        INPUT_BAM = input.bam.vec,
        EXPERIMENT_NAME = paste0( "metpeak_",group.id )
)
bed_name <- paste0( "metpeak_",group.id ,"/peak.xls")
output_bed_name <- paste0("metpeak_",group.id,"_normalized.bed") #peak.bed
bed12.to.bed6 <- paste0("awk 'BEGIN{OFS=\"\t\"}NR>1{print $1,$2,$3,$1\":\"$2\"-\"$3,$5,$6,$7,$8,$9,$10,$11,$12}' ", bed_name," | bed12ToBed6 -i | awk 'BEGIN{FS=\"\t\";OFS=\"\t\"}{print $1,$2,$3,$4,$5}'> ", output_bed_name)
system(bed12.to.bed6)