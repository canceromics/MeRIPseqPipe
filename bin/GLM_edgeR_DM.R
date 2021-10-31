#!/bin/Rscript
## Rscript GLM_edgeR_DM.R <desginfile> <compare_str>
### designfile: Sample_id, Input_filename, IP_filename, group_id
### compare_str: Compairision design (eg: A_vs_B)

#####GLM model###
#edgeR
library(edgeR)

#load data
args <- commandArgs(T)
designfile <- args[1]
compare_str <- args[2]
annotation.file <- args[3]
input.count.matrix.file <- args[4]

designtable <- read.csv(designfile, head = TRUE, stringsAsFactors=FALSE, colClasses = c("character"), check.names=F)
design.matrix <- as.data.frame(designtable$Group)
rownames(design.matrix) <- designtable$Sample_ID
colnames(design.matrix) <- "Condition"

# Get the information of groups from compare_str
if(length(unique(design.matrix$Condition)) < 2){
  stop( "The count of Group is less than two, please check your designfile.")
}else if( compare_str == "two_group" ){
  # Get the information without compare_str beacause of only two groups
  group_id_1 <- unique(design.matrix$Condition)[1]
  group_id_2 <- unique(design.matrix$Condition)[2]
}else{
  # Running MeTDiff quantification with compare_str
  group_id_1 <- strsplit(as.character(compare_str), "_vs_")[[1]][1]
  group_id_2 <- strsplit(as.character(compare_str), "_vs_")[[1]][2]
}
design.matrix <- subset(design.matrix, Condition == group_id_1 | Condition == group_id_2 )
design.matrix$Condition <- factor(design.matrix$Condition, levels = c(group_id_1,group_id_2), labels = c("control","treatment"))
filelist = list.files(path = ".",pattern = ".count",full.names = T)
## Generate the matrix of peaks count
peaks.count.list <- NULL
for(sample_id in row.names(design.matrix)){
  input_count_file <- grep(paste0("[.]",sample_id,"[.]input"),filelist,value = TRUE)
  input_count_table <- read.table(file = input_count_file, sep = "\t", row.names = NULL,header = T)

  ip_count_file <- grep(paste0("[.]",sample_id,"[.]ip"),filelist,value = TRUE)
  ip_count_table <- read.table(file = ip_count_file, sep = "\t", row.names = NULL, header = T)
  rpkm <- cbind(input_count_table[,5],ip_count_table[,5])
  colnames(rpkm) <- c(paste0(sample_id,".input"),paste0(sample_id,".ip"))
  peaks.count.list <- cbind(peaks.count.list,rpkm)
}
rownames(peaks.count.list) <- ip_count_table$PeakName

## generate design matrix
design.matrix$m6A <- "input"
design.matrix$sample_id <- paste0(rownames(design.matrix),".input")
design.matrix_ip <- design.matrix
design.matrix_ip$m6A <- "IP"
design.matrix_ip$sample_id <- paste0(rownames(design.matrix_ip),".ip")
design.matrix <- rbind(design.matrix,design.matrix_ip)
rownames(design.matrix) <- design.matrix$sample_id
design.matrix$m6A <- factor(design.matrix$m6A)
design.matrix <- design.matrix[colnames(peaks.count.list),]

run.edger <- function(cnts,meta){
  #add count filter?
  er.design <- model.matrix(~meta$Condition+meta$m6A+meta$Condition*meta$m6A)
  er.dgelist <- edgeR::DGEList(counts=cnts,group=meta$Condition) 
  er.dgelist <- edgeR::estimateDisp(er.dgelist, design=er.design)
  er.fit <- edgeR::glmFit(er.dgelist, er.design)
  er.lrt <- edgeR::glmLRT(er.fit, coef=4)
  #hist(er.lrt$table$PValue) er.lrt$table$logFC,
  results <- er.lrt$table
  results$padj <- p.adjust(results$PValue,"BH")
  colnames(results) <- c("log2FC","logCPM","LR","pvalue","padj")
  return(results)
}
run.deseq2.4l2fc <- function(cnts,meta,label){
  dds <- DESeq2::DESeqDataSetFromMatrix(cnts,meta,formula(~Condition))
  dds$Condition <- factor(dds$Condition, levels=c('control','treatment'))
  gene.col2check <- meta$Condition
  dds$Condition <- droplevels(dds$Condition)
  gene.deseq <- DESeq2::DESeq(dds)
  gene.deseq <- DESeq2::results(gene.deseq)
  gene.results <- gene.deseq[,c("log2FoldChange","pvalue","padj")]
  colnames(gene.results) <- paste0(label,c(".l2fc",".p",".padj"))
  return(gene.results)
}

results <- run.edger(peaks.count.list,design.matrix)
ip.peaks.count <- read.table(annotation.file,sep = "\t", row.names = 1, check.names = F)
input.gene.count <- read.table(input.count.matrix.file,sep = " ",row.names = 1, check.names = F)
colnames(input.gene.count) <- paste0(colnames(input.gene.count),".input")
peaks.de <- run.deseq2.4l2fc(peaks.count.list[,rownames(design.matrix)[design.matrix$m6A == "IP"]],
                             design.matrix[rownames(design.matrix)[design.matrix$m6A == "IP"],],'peak')
gene.de <- run.deseq2.4l2fc(input.gene.count[,rownames(design.matrix)[design.matrix$m6A == "input"]],
                            design.matrix[rownames(design.matrix)[design.matrix$m6A == "input"],],'gene')
peaks.de$gene_id <- ip.peaks.count[rownames(peaks.de),"ID"]
peaks.de$gene.l2fc <- gene.de[peaks.de$gene_id,]$gene.l2fc
peaks.de$diff.l2fc <- peaks.de$peak.l2fc - peaks.de$gene.l2fc
results$diff.l2fc <-  peaks.de[rownames(results),]$diff.l2fc
write.table(results,file = paste0("edgeR_diffm6A_",group_id_1, "_",group_id_2,".txt") ,sep = "\t",quote = F)

