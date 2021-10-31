#!/bin/Rscript
## Rscript edgeR.R <designfile> <compare_str> eg. Rscript edgeR.R designfile_single.txt T_vs_N
### designfile: Sample_id, Input_filename, IP_filename, group_id
### compare_str: Compairision design (eg: A_vs_B)
library("edgeR")
args<-commandArgs(T) 
designfile <- args[1]
compare_str <- as.character(args[2])

designtable <- read.csv(designfile,head = TRUE,stringsAsFactors=FALSE, colClasses = c("character"))
# Running edgeR by compare.file
## while there are only 2 groups, running edgeR without compare.file
if(length(unique(designtable$Group)) < 2){
  stop( "The count of Group is less than two, please check your designfile.")
}else if( compare_str == "two_group" ){
  # Running edgeR without compare_str beacause of only two groups
  ## Combine expression matrix
  group_id_1 <- unique(designtable$Group)[1]
  group_id_2 <- unique(designtable$Group)[2]
}else{
  # Running edgeR with compare_str 
  ## Combine expression matrix
  group_id_1 <- strsplit(as.character(compare_str), "_vs_")[[1]][1]
  group_id_2 <- strsplit(as.character(compare_str), "_vs_")[[1]][2]
}
control_database = read.table(paste0("htseq_group_", group_id_1, "_input.count"), header = TRUE, row.names = 1)
treated_database = read.table(paste0("htseq_group_", group_id_2, "_input.count"), header = TRUE, row.names = 1)
combined_database <- cbind(control_database,treated_database)
group <- factor(c(rep(group_id_1,ncol(control_database)), rep(group_id_2,ncol(treated_database)))) #setting factors
y <- DGEList(counts=combined_database,group=group)
rownames(y) <- rownames(combined_database)
y <- calcNormFactors(y)
design <- model.matrix(~group)
y <- estimateDisp(y,design)
#To perform likelihood ratio tests:
fit <- glmFit(y,design)
lrt <- glmLRT(fit,coef=2)
topTags(lrt)
### set output_name
lrt$table$padj <- p.adjust(lrt$table$PValue,"BH")
lrt.res <- lrt$table[order(lrt$table$padj),]
colnames(lrt.res) <- c("log2FoldChange","logCPM","LR","pvalue","padj")
output_name <- paste0("edgeR_group_",group_id_1, "_",group_id_2)
write.csv(combined_database, file = paste0(output_name,".matirx") )
write.csv(lrt.res, file = paste0(output_name, ".csv"))
