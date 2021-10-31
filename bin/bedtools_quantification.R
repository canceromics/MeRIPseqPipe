#!/bin/Rscript
## Rscript bedtools_quantification.R <designfile> <bam_stat_summary_file>
### the content of bam_stat_summary_file: example.bam TOTAL_READS
### designfile: Sample_id, Input_filename, IP_filename, group_id
args <- commandArgs(T)
designfile <- args[1]
bam_stat_summary <- args[2]
designtable <- read.csv(designfile,head = TRUE,stringsAsFactors=FALSE, colClasses = c("character"))
bam_stat_table <- read.table(bam_stat_summary,row.names = 1)
filelist =list.files(path = "./",pattern = ".count")
## Generate the quantificative value of peaks referred to RPKM
rpkm_peaks_list <- NULL
#rpkm_peaks_list1 <- NULL
#rpkm_peaks_list2 <- NULL
for(sample_id in designtable$Sample_ID){
  input_count_file <- grep(paste0("[.]",sample_id,"[.]input"),filelist,value = TRUE)
  input_count_table <- read.table(file = input_count_file, sep = "\t", row.names = NULL,header = T)
  bam_stat_index = grep(paste0("^",sample_id,"[.]input"),rownames(bam_stat_table))
  input_rpkm =  apply(input_count_table,1,function(x) (as.numeric(x[5])/(as.numeric(x[3])-as.numeric(x[2]))*1000/bam_stat_table[bam_stat_index,]*1000000))
  
  ip_count_file <- grep(paste0("[.]",sample_id,"[.]ip"),filelist,value = TRUE)
  ip_count_table <- read.table(file = ip_count_file, sep = "\t", row.names = NULL, header = T)
  bam_stat_index = grep(paste0("^",sample_id,"[.]ip"),rownames(bam_stat_table))
  ip_rpkm =  apply(ip_count_table,1,function(x) (as.numeric(x[5])/(as.numeric(x[3])-as.numeric(x[2]))*1000/bam_stat_table[bam_stat_index,]*1000000))
  
  rpkm <- as.matrix((ip_rpkm+1)/(input_rpkm+1))
  #rpkm1<- as.matrix((ip_rpkm)/(input_rpkm+1))
  #rpkm2<- as.matrix((ip_rpkm)/(input_rpkm+ip_rpkm))
  colnames(rpkm)[1] <- sample_id
  #colnames(rpkm1)[1] <- sample_id
  #colnames(rpkm2)[1] <- sample_id
  rpkm_peaks_list <- cbind(rpkm_peaks_list,rpkm)
  #rpkm_peaks_list1 <- cbind(rpkm_peaks_list1,rpkm1)
  #rpkm_peaks_list2 <- cbind(rpkm_peaks_list2,rpkm2)	
}
rownames(rpkm_peaks_list) <- input_count_table[,4] 
#rownames(rpkm_peaks_list1) <- input_count_table[,4] 
#rownames(rpkm_peaks_list2) <- input_count_table[,4] 
write.table(rpkm_peaks_list,sep = "\t",file = "bedtools_quantification.matrix",quote = F)
#write.table(rpkm_peaks_list1,sep = "\t",file = "bedtools_quantification.inputadd.matrix",quote = F)
#write.table(rpkm_peaks_list2,sep = "\t",file = "bedtools_quantification.alladd.matrix",quote = F)
