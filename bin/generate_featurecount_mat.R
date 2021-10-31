#!/bin/Rscript
## Rscript get_htseq_matrix.R designfile THREAD_NUM eg. Rscript get_htseq_matrix.R designfile_single.txt 10
## designfile: filename, control_or_treated, input_or_ip, group(default 0 is CONTROL_SITUATION else are TREATED_SITUATION)

library(parallel)
library(data.table)
args<-commandArgs(T)
designfile <- args[1]
THREAD_NUM <- as.numeric(args[2])

designtable <- read.csv(designfile,header = TRUE,stringsAsFactors=FALSE, colClasses = c("character"))
#Generate gene count matrix
fc.files <- list.files("./",pattern = ".txt$")
mclapply(unique(designtable$Group),function(x){
  group_id <- x
  group.input.count.mat <- NULL
  for(pc in grep(paste0(".input_",group_id,"[.]bam"),fc.files,value = TRUE)){
    pc.exp <- fread(pc,sep = "\t")[,c(1,7)]
    if(is.null(group.input.count.mat)){
      group.input.count.mat <- pc.exp
    }else{
      group.input.count.mat <- merge(group.input.count.mat,pc.exp,by = c("Geneid"))
    }
  }
  #parsing samplenames
  output_pattern = paste0("htseq_group_",group_id)  #添加aligner
  fwrite(group.input.count.mat, file = paste0(output_pattern,"_input.count"), sep = "\t")
},
mc.cores = THREAD_NUM
)
group.mat.list = grep("htseq",list.files(path = "./",pattern = "input.count"), value = T)
expression.matrix <- NULL
for( file in group.mat.list ){
  tmp.expression.table <- as.matrix(read.table(file, header = TRUE, row.names = 1, check.names=F))
  expression.matrix <- cbind(expression.matrix, tmp.expression.table)
}
colnames(expression.matrix) <- as.matrix(lapply(strsplit(colnames(expression.matrix),".input"), function(x){ x[1]}))
write.table(expression.matrix,file = "expression.matrix",quote=F)
