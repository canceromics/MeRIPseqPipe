#!/bin/Rscript
## Rscript QC_Peaks_Report.R <desginfile> <peakMerged_mode>
### designfile: Sample_id, Input_filename, IP_filename, group_id
library(ggplot2)
library(ggseqlogo)
library(reshape2)
options(stringsAsFactors = FALSE)

args <- commandArgs(T) 
#args <- c("formatted_designfile.txt", "mspc", "group","QCPeaksPlot.RData")
designfile <- args[1] #"formatted_designfile.txt"
peakMerged.mode <- args[2]#rank
peakCalling.mode <- args[3]#"group"
output.Rdata <- args[4]#"QCPeaksPlot.RData"
designtable <- read.csv(designfile, head = TRUE, colClasses = c("character"))

## Peaks Distribution
pdf(file = paste0("distribution.plot_",peakMerged.mode,".pdf"),paper = "USr")
PlotPeaksDitr <- function(files.list, suffix = "[.]anno[.]txt"){
  distribute_df <- NULL
  for( file in files.list ){
    anno.table <- read.table(file, header=F, sep="\t", quote="", stringsAsFactors = F)[,c(1,2,3,15,11,12,13,14,17)]
    colnames(anno.table) <- c("Chr","ChrStart","ChrEnd","ID","Gene_symbol","Coding","Location","Relative_distance","RNA_type")
    peak.freq = c(as.numeric(as.vector(anno.table[which(anno.table[,7]=="5UTR"),8])),
                  as.numeric(as.vector(anno.table[which(anno.table[,7]=="CDS"),8]))+100,
                  as.numeric(as.vector(anno.table[which(anno.table[,7]=="3UTR"),8]))+200)
    freq = data.frame(Freq = peak.freq, group = strsplit(file,suffix)[[1]])
    distribute_df = rbind(distribute_df, freq)
  }
  ggplot(distribute_df, aes(x=Freq, colour = group))+
    geom_line(stat = "density", size=1, adjust = 0.8)+
    scale_x_continuous(breaks = c(50,150,250), labels = c("5'UTR", "CDS", "3'UTR"))+ #axis labels
    labs(y="m6A coding peak density",x="Region of gene")+
    geom_vline(xintercept = c(100,200), linetype = "dashed")+
    theme_bw()+
    theme(panel.grid =element_blank(),#remove grid line
          axis.title.x = element_text(size = 20, angle = 0, face = "plain", colour = "black"),
          axis.title.y = element_text(size = 20, angle = 90, face = "plain", colour = "black"),
          axis.text.x = element_text(size = 15,colour = "black"),
          axis.text.y = element_text(size = 15,colour = "black"),
          aspect.ratio=1,
          axis.ticks.x = element_blank()) #remove ticks
}

anno.files.list <- dir(pattern = "[.]anno[.]txt")
### barplot
total.distribute <- NULL
for( file in anno.files.list ){
  anno.table <- read.table(file, header=F, sep="\t", quote="", stringsAsFactors = F)[,c(1,2,3,15,11,12,13,14,17)]
  colnames(anno.table) <- c("Chr","ChrStart","ChrEnd","ID","Gene_symbol","Coding","Location","Relative_distance","RNA_type")
  anno.table[anno.table$Location == "CDS" & anno.table$Relative_distance >= 95,7] <- "Stop Codon"
  anno.table[anno.table$Location == "3UTR" & anno.table$Relative_distance <= 5,7] <- "Stop Codon"
  group.name <- strsplit(file,"[.]anno[.]txt")[[1]]
  freq = data.frame(Location = anno.table[anno.table$Coding == "coding",7], group = group.name)
  total.distribute = rbind(total.distribute, freq)
}
distribute.table <- melt(table(total.distribute))
distribute.table$Location <- factor(distribute.table$Location,levels = c("intron","3UTR","Stop Codon","CDS","5UTR"))
col <- c('plum1','pink2','#58B2DC',"#51A8DD","#005CAF")
distribute.barplot <- ggplot(distribute.table,aes(group, value, fill = Location)) +
                        geom_bar(stat="identity",position = 'fill') + coord_flip() +
                        ggtitle("Peaks Distribution") +
                        scale_y_continuous(expand = c(0, 0)) +
                        guides(fill = guide_legend(reverse = TRUE)) +
                        scale_fill_brewer() +
                        theme(panel.grid =element_blank(), #remove grid line
                              title = element_text(size = 15, angle = 0, face = "plain", colour = "black"),
                              axis.text.x = element_text(size = 12,colour = "black"),
                              axis.text.y = element_text(size = 12,colour = "black"),
                              panel.background = element_rect(fill = "transparent",colour = NA),
                              axis.title = element_blank(),
                              axis.ticks.x = element_blank()) #remove ticks
distribute.barplot.count <- ggplot(distribute.table,aes(group, value, fill = Location)) +
                            geom_bar(stat="identity") + coord_flip() +
                            ggtitle("Peaks Distribution") +
                            scale_y_continuous(expand = c(0, 0)) +
                            guides(fill = guide_legend(reverse = TRUE)) +
                            scale_fill_brewer() +
                            theme(panel.grid =element_blank(), #remove grid line
                                  title = element_text(size = 15, angle = 0, face = "plain", colour = "black"),
                                  axis.text.x = element_text(size = 12,colour = "black"),
                                  axis.text.y = element_text(size = 12,colour = "black"),
                                  panel.background = element_rect(fill = "transparent",colour = NA),
                                  axis.title = element_blank(),
                                  axis.ticks.x = element_blank()) #remove ticks
print(distribute.barplot)
print(distribute.barplot.count)
### Curve
sample.plots.list <- NULL
sample.list <- if(peakCalling.mode == "group") designtable$Group else designtable$Sample_ID
for( sample in sample.list ){
  sample.files.list <- grep(paste0("_",sample,"_normalized"), anno.files.list, value = T)
  sample.plots.list[[sample]] <- PlotPeaksDitr(sample.files.list, "_normalized[.]anno[.]txt")
  print(sample.plots.list[[sample]])
}
merged.files.list <- grep("merged", anno.files.list, value = T)
merged.plot <- PlotPeaksDitr(merged.files.list)
print(merged.plot)
dev.off()

## Peaks' motif
pdf(file = paste0("motif.plot_",peakMerged.mode,".pdf"),paper = "USr")
ggplot2.multiplot <- function(..., plotlist=NULL, cols=2) {
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  numPlots = length(plots)
  
  # Make the panel
  plotCols = cols                          # Number of columns of plots
  plotRows = ceiling(numPlots/plotCols) # Number of rows needed, calculated from # of cols
  
  # Set up the page
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(layout = grid::grid.layout(plotRows, plotCols)))
  vplayout <- function(x, y)
    grid::viewport(layout.pos.row = x, layout.pos.col = y)
  
  # Make each plot, in the correct location
  for (i in 1:numPlots) {
    curRow = ceiling(i/plotCols)
    curCol = (i-1) %% plotCols + 1
    print(plots[[i]], vp = vplayout(curRow, curCol))
  }
  
}
motif_plot <- function(motif, pval, rank){
  ggplot()+
    geom_logo(motif, method = "probability")+
    annotate("text", x=ncol(motif)-0.5, y=1.5, label=paste0("p = ",pval),size = 5)+
    ggtitle(rank)+
    theme(plot.title = element_text(hjust = 0, size = 6))+
    theme_logo()
}
QC.motif.filelist = dir(".",pattern = "motif[1,2,3].motif",recursive = TRUE)
QC.motif.list <- NULL
QC.motif.pvalue <- NULL
motif.peakfiles <- unique(unlist(lapply(strsplit(QC.motif.filelist,"_homer/homerResults"), function(x){x[1]})))
for( peakfile.name in motif.peakfiles ){
  group.motif.list <- NULL
  for (file in grep(paste0(peakfile.name, "_homer/homerResults"), QC.motif.filelist, value = T) ){
    motif_matrix <- read.delim(file,header = F,sep = "\t", check.names=F)
    motif_pvalue <- strsplit(motif_matrix[1,6], split = ":")[[1]][4]
    motif_matrix <- motif_matrix[-1,c(-5,-6)]
    colnames(motif_matrix) <- c("A","C","G","T")
    rownames(motif_matrix) <- c(1:nrow(motif_matrix))
    motif_matrix <- as.matrix(t(motif_matrix))
    motif_name <- strsplit(strsplit(file,split = c("_homer/homerResults/"))[[1]][2],split = "[.]motif")[[1]][1]
    group.motif.list[[motif_name]] <- motif_plot(t(apply(motif_matrix, 1, function(x)as.numeric(x))), motif_pvalue, paste0(peakfile.name,"_",motif_name))
  }
  QC.motif.list[[peakfile.name]] <- group.motif.list
  ggplot2.multiplot(plotlist = QC.motif.list[[peakfile.name]] ,cols = 1)
}
dev.off()
save(distribute.barplot.count,distribute.barplot,sample.plots.list,merged.plot,QC.motif.list,file = output.Rdata)
