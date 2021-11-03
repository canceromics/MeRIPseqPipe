#!/bin/Rscript
## Rscript get_featurecount_matrix.R designfile THREAD_NUM gtf eg. Rscript get_featurecount_matrix.R designfile_single.txt 10 gtf_file
## designfile: filename, control_or_treated, input_or_ip, group(default 0 is CONTROL_SITUATION else are TREATED_SITUATION)
library(GenomicFeatures)
library(parallel)
library(data.table)
args<-commandArgs(T)
designfile <- args[1]
THREAD_NUM <- as.numeric(args[2])
gtf <- args[3]
aligner <- args[4]

designtable <- read.csv(designfile,header = TRUE,stringsAsFactors=FALSE, colClasses = c("character"))
#Generate gene count matrix
txdb <- makeTxDbFromGFF(gtf,format = 'gtf')
exonsgene <- exonsBy(txdb,by = 'gene')
exonsgene_len <- lapply(exonsgene,function(x){sum(width(reduce(x)))})
len <- t(as.data.frame(exonsgene_len))

count2fpkm <- function(counts, effLen)
{
  N <- sum(counts)
  exp( log(counts) + log(1e9) - log(effLen) - log(N) )
}

count2tpm <- function(counts, effLen)
{
  rate <- log(counts) - log(effLen)
  denom <- log(sum(exp(rate)))
  exp(rate - denom + log(1e6))
}

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
  output_pattern = paste0("featurecount_",aligner,"_group_",group_id)  #添加aligner
  fwrite(group.input.count.mat, file = paste0(output_pattern,"_input.count"), sep = "\t")
},
mc.cores = THREAD_NUM
)
group.mat.list = grep("featurecount",list.files(path = "./",pattern = "input.count"), value = T)
expression.matrix <- NULL
for( file in group.mat.list ){
  tmp.expression.table <- as.matrix(read.table(file, header = TRUE, row.names = 1, check.names=F))
  expression.matrix <- cbind(expression.matrix, tmp.expression.table)
}
colnames(expression.matrix) <- as.matrix(lapply(strsplit(colnames(expression.matrix),".input"), function(x){ x[1]}))
write.table(expression.matrix,file = paste0("expression.",aligner,".count.matrix"),quote=F)

expression.matrix <- as.data.frame(expression.matrix)
fpkm_matrix <- count2fpkm(expression.matrix[,1],len)
for (i in 2:length(expression.matrix)) {
  temp <- count2fpkm(expression.matrix[,i],len)
  fpkm_matrix <- data.frame(fpkm_matrix,temp)
}
colnames(fpkm_matrix) <- colnames(expression.matrix)
write.table(fpkm_matrix,file = paste0("expression.",aligner,".fpkm.matrix"),quote = F)

tpm_matrix <- count2fpkm(expression.matrix[,1],len)
for (i in 2:length(expression.matrix)) {
  temp <- count2tpm(expression.matrix[,i],len)
  tpm_matrix <- data.frame(tpm_matrix,temp)
}
colnames(tpm_matrix) <- colnames(expression.matrix)
write.table(tpm_matrix,file = paste0("expression.",aligner,".tpm.matrix"),quote = F)