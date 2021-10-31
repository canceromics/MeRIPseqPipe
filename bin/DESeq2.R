#!/bin/Rscript
## Rscript DESeq2.R <designfile> <compare_str>  eg. Rscript DESeq2.R designfile_single.txt T_vs_N 
### designfile: Sample_id, Input_filename, IP_filename, group_id
### compare_str: Compairision design (eg: A_vs_B)

library(DESeq2)
args<-commandArgs(T) 
designfile <- args[1]
compare_str <- as.character(args[2])

designtable <- read.csv(designfile,head = TRUE,stringsAsFactors=FALSE, colClasses = c("character"))
# Running DEseq2 by compare_str
## while there are only 2 groups, running DEseq2 without compare_str
if( length(unique(designtable$Group)) < 2 ){
  stop( "The count of Group is less than two, please check your designfile.")
}else if( compare_str == "two_group" ){
  # Running DESeq2 without compare_str beacause of only two groups
  ## Combine expression matrix
  group_id_1 <- unique(designtable$Group)[1]
  group_id_2 <- unique(designtable$Group)[2]
}else{
  # Running DESeq2 with compare_str
  ## Combine expression matrix
  group_id_1 <- strsplit(as.character(compare_str), "_vs_")[[1]][1]
  group_id_2 <- strsplit(as.character(compare_str), "_vs_")[[1]][2]
}
control_database = read.table(paste0("htseq_group_", group_id_1, "_input.count"), header = TRUE, row.names = 1, check.names = FALSE)
treated_database = read.table(paste0("htseq_group_", group_id_2, "_input.count"), header = TRUE, row.names = 1, check.names = FALSE)
combined_database <- cbind(control_database,treated_database)
condition <- factor(c(rep(group_id_1,ncol(control_database)), rep(group_id_2,ncol(treated_database)))) #setting factors
### assign gene names
colData <- data.frame(row.names=colnames(combined_database), group = condition)
dds <- DESeqDataSetFromMatrix(countData = combined_database,colData = colData,design = ~ group)
rownames(dds) <- rownames(combined_database)
#dds <- dds[ rowSums(counts(dds)) > 1, ]
dds <- DESeq(dds)
## FoldChange = group_id_2 / group_id_1
res <- results(object = dds, contrast = c("group",group_id_2,group_id_1))
table(res$padj <0.05)
res <- res[order(res$padj),]
#resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)),by="row.names",sort=FALSE)
#resdata2=resdata[resdata$log2FoldChange > 1|resdata$log2FoldChange < -1, ]
### set output_name
output_name <- paste0("DESeq2_group_",group_id_1, "_",group_id_2)
write.csv(res, file = paste0(output_name, ".csv"))
#write.csv(resdata2,file = paste0(output_name, "_log2.csv"),row.names =FALSE)
