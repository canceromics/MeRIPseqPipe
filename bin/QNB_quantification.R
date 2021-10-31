# Rscript QNB_quantification.R designfile
library("QNB")
args <- commandArgs(T)
designfile <- args[1]
## read designfile to get the name of samples
designtable <- read.csv(designfile,head = TRUE,stringsAsFactors=FALSE, colClasses = c("character"))
filelist = list.files(path = "./",pattern = ".count")
rpkm_peaks_list <- NULL
for(sample_id in designtable$Sample_ID){
  ## generate the dataframe of peak count
  input.count <- c()
  input.names <- c()
  input.samples <- c()
  for(input in grep(sample_id, grep("[.]input",filelist,value = TRUE), value = T)){
    input.exp <- read.table(input,header=T,sep="\t",row.names= NULL,quote = "")
    input.count <- cbind(input.count,input.exp[,5])
    input.names <- input.exp[,4] #peaks name
    input.samples <- c(input.samples,input) #samples name
  }
  colnames(input.count) <- input.samples
  rownames(input.count) <- input.names
  ip.count <- c()
  ip.names <- c()
  ip.samples <- c()
  for(ip in grep(sample_id, grep("[.]ip",filelist,value = TRUE), value = T)){
    ip.exp <- read.table(ip,header=T,sep="\t",row.names= NULL,quote = "")
    ip.count <- cbind(ip.count,ip.exp[,5])
    ip.names <- ip.exp[,4] #peaks name
    ip.samples <- c(ip.samples,ip) #samples name
  }
  colnames(ip.count) <- ip.samples
  rownames(ip.count) <- ip.names
  ## Run the QNB to generate quantificative value per sample
  result <- qnbtest(ip.count, ip.count, input.count, input.count, mode="blind")
  sample_quantification <- as.matrix(result$p.treated)
  colnames(sample_quantification) <- sample_id
  rpkm_peaks_list <- cbind(rpkm_peaks_list,sample_quantification)
  rownames(rpkm_peaks_list) <- rownames(ip.count)
}
write.table(rpkm_peaks_list,sep = "\t",file = "QNB_quantification.matrix",quote = F)
