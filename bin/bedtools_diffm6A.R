#!/bin/Rscript
## Rscript bedtools_diffm6A.R <desginfile> <quantification_matrix_file> <compare_str>
### designfile: Sample_id, Input_filename, IP_filename, group_id
### compare_str: Compairision design (eg: A_vs_B)
args <- commandArgs(T)
designfile <- args[1]
quantification_matrix_file <- args[2]
compare_str <- as.character(args[3])

designtable <- read.csv(designfile, head = TRUE, stringsAsFactors=FALSE, colClasses = c("character"))
quantification_matrix = read.table(quantification_matrix_file ,sep = "\t",header = T, row.names = 1)
# generate design matrix
design.matrix <- as.data.frame(designtable$Group)
rownames(design.matrix) <- designtable$Sample_ID
colnames(design.matrix) <- "Type"

# Wilcoxon test for the vector of two groups
row_wilcox <- function(design.matrix,group_id_1,group_id_2,x,test_mode=""){
  group1 <- as.character(rownames(subset(design.matrix,Type==group_id_1))) 
  group2 <- as.character(rownames(subset(design.matrix,Type==group_id_2)))
  if (test_mode=="paired"){
    res_wix0 <- wilcox.test(x[which(rownames(design.matrix)%in%group1)],x[which(rownames(design.matrix)%in%group2)], paired = T)
  } else {
    res_wix0 <- wilcox.test(x[which(rownames(design.matrix)%in%group1)],x[which(rownames(design.matrix)%in%group2)])
  }
  res_wix0$log2FC = log2(mean(x[which(rownames(design.matrix)%in%group2)])/mean(x[which(rownames(design.matrix)%in%group1)]))
  res_wix <- c(log2FC=res_wix0$log2FC,pvalue=res_wix0$p.value,statistic=res_wix0$statistic) 
  return(res_wix)
}

# Get the information of groups from compare_str
if(length(unique(design.matrix$Type)) < 2){
  stop( "The count of Group is less than two, please check your designfile.")
}else if( compare_str == "two_group" ){
  # Get the information without compare_str beacause of only two groups
  group_id_1 <- unique(design.matrix$Type)[1]
  group_id_2 <- unique(design.matrix$Type)[2]
}else{
  # Running MeTDiff quantification with compare_str
  group_id_1 <- strsplit(as.character(compare_str), "_vs_")[[1]][1]
  group_id_2 <- strsplit(as.character(compare_str), "_vs_")[[1]][2]
}
cat("peak number ",dim(quantification_matrix)[1],"\n")

# Run the Wilcoxon test for the quantifacative value of every peak
test_mode=""
res_wix_lst <- apply(na.omit(quantification_matrix[,3:ncol(quantification_matrix)]),1,function(x){row_wilcox(design.matrix,group_id_1,group_id_2,x,test_mode)})
res_wix_lst = as.data.frame(t(res_wix_lst)) 
res_wix_lst$padj = p.adjust(res_wix_lst$pvalue,method = "BH")
res_wix_lst$BY = p.adjust(res_wix_lst$pvalue,method = "bonferroni")
cat("DM peaks pvalue(0.05)",sum(res_wix_lst$pvalue <=0.05),"\n")
cat("DM peaks FDR(0.05)",sum(res_wix_lst$padj <=0.05),"\n")
output_name <- paste0("bedtools_diffm6A_",group_id_1, "_", group_id_2)
write.table(res_wix_lst, file = paste0(output_name,".txt"), sep = "\t", quote = F)

