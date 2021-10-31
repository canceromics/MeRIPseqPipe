# required for merge_peaks_by_rank.sh
library(RobustRankAggreg)
args<-commandArgs(T)
bedlist <- read.table(args[1],header = F,sep = "\t",stringsAsFactors = F, na.strings = "")
len_of_bed <- as.numeric(args[2])
out_name <- as.character(args[3])
bedlist2 <- as.list(NULL)
for (i in c(1:ncol(bedlist))){
  if (TRUE %in% is.na(bedlist[,i])){
    sub <- which(is.na(bedlist[,i]))
    bedlist2[[i]] <- bedlist[-sub,i]
  }
  else{
    bedlist2[[i]] <-bedlist[,i]
  }
}
mergepeak <- aggregateRanks(glist = bedlist2, N = len_of_bed)
sub <- which(as.numeric(mergepeak$Score)==1)
mergepeak <- mergepeak[-sub,]
merged.bed <- apply(mergepeak, 1 ,function(x){
  peak.info <- as.vector(as.matrix(x))
  peak.region = unlist(strsplit(strsplit(as.character(peak.info[1]),split = ":" )[[1]],split = "-"))
  x = c(peak.region,as.character(peak.info[1]),peak.info[2])
})
## 
merged.bed <- t(merged.bed)
write.table(merged.bed,file = out_name,sep = "\t",quote = FALSE,row.names = FALSE,col.names = FALSE)
