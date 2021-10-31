#!/bin/Rscript
## Rscript arranged_result.R <designfile> <comparefile> <diffm6A_mode> <expression_mode> <mergePeaks_mode>
### designfile: Sample_id, Input_filename, IP_filename, group_id
### compare_str: Compairision design (eg: A_vs_B)
args <- commandArgs(T)
#args <- c("macs2_MeTDiff_DESeq2_arranged_results_2019-12-11.m6APipe", "DiffReport.RData")
m6APipe.data <- args[1]#"formatted_designfile.txt"
output.Rdata <- args[2]#"compare_info"

library(pheatmap)
library(ggplot2)
library(ggrepel)
library(grid)
library(reshape2)
load(m6APipe.data)

draw_colnames_90 <- function (coln, ...) {
  m = length(coln)
  x = (1:m)/m - 1/2/m
  grid.text(coln, x = x, y = unit(0.96, "npc"), vjust = .5, 
            hjust = 1, rot = 90, gp = gpar(...)) ## 注意缺省值为 'hjust=0' 和'rot=270'
}
assignInNamespace(x="draw_colnames", value="draw_colnames_90",
                  ns=asNamespace("pheatmap"))

heatmap_dm <- function(mat,coldt){
  pheatmap(mat, cluster_rows=FALSE, show_rownames=F, cluster_cols=FALSE, annotation_col=coldt, 
           main = "Heatmap of Different Methylation", scale = "row")
}

heatmap_de <- function(mat, coldt){
  pheatmap(mat, cluster_rows=FALSE, show_rownames=F, cluster_cols=FALSE, annotation_col=coldt, 
           color = colorRampPalette(c(rep('#1C2B6F',1),'black', rep('#E31E26',1)))(50),
           main = "Heatmap of Different Expression", scale = "row")
}

ECDF_plot <- function(df,value_var,group,plot_title="",test_result=""){
  if (test_result !=""){
    if (test_result$p.value <=0){
      test_anno = paste(test_result$method,"\n",names(test_result$statistic)," = ",signif(test_result$statistic, 3),
                        "\nP Value < 2.2e-16",sep = "")
    } else {
      test_anno = paste(test_result$method,"\n",names(test_result$statistic)," = ",signif(test_result$statistic, 3),
                        "\nP Value= ",signif(test_result$p.value,3),sep = "")
    }
  }else{
    test_anno = ""
  }
  p <-  ggplot(df,aes(x=value_var,group=group,color=group))+theme_test()+
    stat_ecdf(size = 1)+theme(legend.position=c(0.85,.15))+
    annotate("text",x=-Inf,y=Inf,vjust=1.5,hjust=-.12,label=test_anno)+ 
    scale_y_continuous(expand = c(0,0))+scale_x_continuous(expand = c(0,0),limits = c(0,1.02))+
    labs(title= plot_title, y="Cumulative fraction" , x = "Peaks intensity")+
    theme(plot.title = element_text(size = 15, angle = 0, face = "plain", colour = "black"),
          axis.title.x = element_text(size = 15, angle = 0, face = "plain", colour = "black"),
          axis.title.y = element_text(size = 15, angle = 90, face = "plain", colour = "black"),
          axis.text.x = element_text(size = 15, angle = 0, face = "plain", colour = "black"),
          axis.text.y = element_text(size = 15, angle = 0, face = "plain", colour = "black"))
  return(p)
}

volcano_plot_dm = function(res, Sample_1 = "A", Sample_2 = "B", lfc = 0.58, pval = 0.05, groupname = ""){
  par(mar = c(5, 6, 5, 5))
  tab = data.frame(logFC = res$log2FC, negLogPval = -log10(res$pvalue)) 
  tab$gene_name = rownames(res)
  tab = na.omit(tab)
  tab$threshold <- "C"
  tab$threshold[tab$logFC >= lfc & tab$negLogPval > -log10(pval)] <- "B"
  tab$threshold[tab$logFC <=-lfc & tab$negLogPval > -log10(pval)] <- "A"
  #tab<-tab%>%mutate(threshold = ifelse(logFC >= lfc & negLogPval > -log10(pval) ,"B", ifelse(logFC<=-lfc & negLogPval > -log10(pval), "A", "C")))
  n_up = length(which(tab$threshold=="B"))
  n_down = length(which(tab$threshold=="A"))
  tab_order = tab[order(tab$negLogPval, decreasing = T),]
  ggplot(tab_order, aes(x=logFC, y=negLogPval)) +
    geom_point(aes(colour = threshold)) +
    scale_colour_manual(values = c("A"= "#619cff", "B"="#f8766d",  "C"= "#c8c8c8"),
                        labels=c(paste("Down: ", n_down, sep=""),paste("Up: ", n_up, sep = "") , "No sig"), name = NULL) +
    geom_hline(aes(yintercept=-log10(pval)), linetype="dashed") +
    geom_vline(aes(xintercept=-lfc), linetype="dashed") +
    geom_vline(aes(xintercept=lfc), linetype="dashed") +
    ggtitle(paste("Volcano Plot of Different Methylation in", groupname))+
    xlab(expression(paste(Log[2], " fold change", sep = ""))) +
    ylab(expression(paste(-Log[10], " adjusted P value", sep = ""))) +
    theme_bw() +
    theme(legend.position = 'top',
          plot.title = element_text(hjust = 0.5))
}

quadrant_plot <- function(quadrant.data, lfc = 0.58 , pval = 0.05, groupname = ""){
    quadrant.data$threshold <- "nosig"
    quadrant.data$threshold[quadrant.data$m6A >= lfc & quadrant.data$exp >= lfc &
                            quadrant.data$m6A.p <= pval & quadrant.data$exp.p <= pval ] <- "Hyper-up"
    quadrant.data$threshold[quadrant.data$m6A >= lfc & quadrant.data$exp <= -lfc &
                            quadrant.data$m6A.p <= pval & quadrant.data$exp.p <= pval ] <- "Hyper-down"
    quadrant.data$threshold[quadrant.data$m6A <= -lfc & quadrant.data$exp >= lfc &
                              quadrant.data$m6A.p <= pval & quadrant.data$exp.p <= pval ] <- "Hypo-up"
    quadrant.data$threshold[quadrant.data$m6A <= -lfc & quadrant.data$exp <= -lfc &
                              quadrant.data$m6A.p <= pval & quadrant.data$exp.p <= pval ] <- "Hypo-down"
    quadrant.data.length <- table(quadrant.data$threshold)
    quadrant.data <- na.omit(quadrant.data)
    quadrant.data$threshold <- factor(quadrant.data$threshold,levels = c("Hyper-up","Hyper-down","Hypo-up","Hypo-down","nosig"))
    ggplot()+
      geom_point(data = quadrant.data,
                 aes_string(x= "exp" ,y="m6A", color="threshold"),size = 1)+
      geom_hline(yintercept = lfc,linetype="dashed")+
      geom_hline(yintercept = -lfc,linetype="dashed")+
      geom_vline(xintercept = lfc,linetype="dashed")+ylab("dmlogfc")+
      geom_vline(xintercept = -lfc,linetype="dashed")+
      scale_x_continuous(limits = c(-5,5))+
      scale_y_continuous(limits = c(-5,5))+
      ggtitle(paste("Quadrant Plot between Methylation and Expression in", groupname))+
      scale_colour_manual(values = c("Hyper-up" = "#7DB9DE", "Hyper-down" = "#D75455", 
                                     "Hypo-up" = "#7BA23F", "Hypo-down" = "#A35E47", 
                                     "nosig"= "#c8c8c8" ),
                          labels = c(paste("Hyper-up: ",   quadrant.data.length["Hyper-up"], sep=""),
                                     paste("Hyper-down: ",   quadrant.data.length["Hyper-down"], sep = ""),
                                     paste("Hypo-up: ",   quadrant.data.length["Hypo-up"], sep=""),
                                     paste("Hypo-down: ",   quadrant.data.length["Hypo-down"], sep = ""),
                                     "No sig"), name = NULL) +
      theme_classic() +
      theme(plot.title = element_text(hjust = 0.5))+
      geom_point(data = quadrant.data[quadrant.data$threshold!="nosig",],
                 aes_string(x= "exp" ,y="m6A", color="threshold"),size = 1.5)

}
matrixcluster <- function(matrixData, cluster_rows = TRUE, cluster_cols = TRUE, cmethod = "complete"){
  if(cluster_rows == TRUE){
    ht <- hclust(dist(matrixData), method = cmethod)          #对行进行聚类
    rowInd <- ht$order                                       #将聚类后行的顺序存为rowInd 
  }else{
    rowInd <- 1:nrow(matrixData)
  }
  
  if(cluster_cols == TRUE){
    ht <- hclust(dist(t(matrixData)), method = cmethod)       #对矩阵进行转置，对原本的列进行聚类
    colInd <- ht$order                                       #将聚类后列的顺序存为colInd 
  }else{
    colInd <- 1:ncol(matrixData)
  }
  
  matrixDataNew <-matrixData[rowInd,colInd]                #将数据按照聚类结果重排行和列
  print(c("ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid"))
  return(matrixDataNew)
}
## dm & de & ecdf
heatmap_dm.list <- NULL
heatmap_de.list <- NULL
volcano_dm.list <- NULL
ecdf.list <- NULL
quadrant.list <- NULL
ecdf.data <- NULL
ecdf.group.data <- melt(m6a.anno.matrix,ID="PeakRegion")
for( group in as.character(compare.list) ){
  group1 = strsplit(group, "_vs_")[[1]][1] 
  group2 = strsplit(group, "_vs_")[[1]][2]
  coldata = subset(as.data.frame(design.matrix), Type==group1|Type==group2)
  coldata$Type = as.factor(coldata$Type)
  ## dm
  dmres = diffm6A.list[[which(names(diffm6A.list)==group)]]
  dmg = subset(dmres, abs(log2FC) > 0.58 & pvalue < 0.05)
  if (nrow(dmg)<=1) dmg = dmres
  matrix.dm = m6a.anno.matrix[,-c(1:3)]
  dm_mat = matrix.dm[dmg$PeakRegion,rownames(coldata)]
  select <- dmg[order(dmg$log2FC, decreasing = TRUE),] 
  dm_mat = log2(dm_mat+1)
  dm_mat = dm_mat[select$PeakRegion,]
  dm_mat = na.omit(dm_mat)
  dm_mat_new <- matrixcluster(dm_mat,cmethod = "single")
  heatmap_dm.list[[group]] <- heatmap_dm(dm_mat_new,coldata)
  ## de
  deres = diffexpression.list[[which(names(diffexpression.list)==group)]]
  deg = subset(deres, abs(log2FoldChange)> 0.58 & pvalue < 0.05)
  if (nrow(deg)<=1) deg = deres
  rownames(deg) = deg$ID
  de_mat = expression.matrix[row.names(deg),rownames(coldata)]
  select <- deg[order(deg$log2FoldChange, decreasing = TRUE), ] 
  de_mat = log2(de_mat+1)
  de_mat = de_mat[rownames(select),]
  de_mat = na.omit(de_mat)
  de_mat_new <- matrixcluster(de_mat,cmethod = "single")
  heatmap_de.list[[group]] <- heatmap_de(de_mat_new,coldata)
  
  ## ecdf
  ecdf.group.data.tmp = subset(ecdf.group.data,variable %in% rownames(coldata))
  ecdf.group.data.tmp$group <- group1
  ecdf.group.data.tmp$group[ecdf.group.data.tmp$variable %in% rownames(coldata)[coldata$Type == group2]] <- group2
  #ecdf.group.data.tmp <- ecdf.group.data.tmp%>%mutate(group = ifelse(variable %in% rownames(coldata)[coldata$Type == group1], group1, group2))
  ecdf.group.data.tmp <- na.omit(ecdf.group.data.tmp)
  ecdf.list[[group]] <- ECDF_plot(ecdf.group.data.tmp,ecdf.group.data.tmp$value,ecdf.group.data.tmp$group)
  ecdf.data <- rbind(ecdf.data,data.frame(data = diffm6A.list[[group]]$log2FC , group = group))
  ## volcano plot
  volcano_dm.list[[group]] <- volcano_plot_dm(diffm6A.list[[group]],groupname = group)
  ## quadrant plot
  diffm6a.results <- diffm6A.list[[group]]
  rownames(diffm6a.results) <- diffm6a.results$PeakRegion
  diffexp.results <- diffexpression.list[[group]]
  rownames(diffexp.results) <- diffexp.results$ID
  quadrant.data <- data.frame(row.names = rownames(diffm6a.results), 
                              m6A = diffm6a.results$log2FC, 
                              m6A.p = diffm6a.results$pvalue,
                              exp = diffexp.results[diffm6a.results$ID,"log2FoldChange"],
                              exp.p = diffexp.results[diffm6a.results$ID,"pvalue"])
  quadrant.list[[group]] <-quadrant_plot(quadrant.data,lfc = 0.58,pval = 0.05,groupname = group)
  
  ## plot
  pdf(file = paste0("heatmap_dm_",group,".pdf"),paper = "USr")
  print(heatmap_dm.list[[group]])
  dev.off()
  pdf(file = paste0("heatmap_de_",group,".pdf"),paper = "USr")
  print(heatmap_de.list[[group]])
  dev.off()
  pdf(file = paste0("volcano_dm_",group,".pdf"),paper = "USr")
  print(volcano_dm.list[[group]])
  dev.off()
  pdf(file = paste0("ecdf_",group,".pdf"),paper = "USr")
  print(ecdf.list[[group]])
  dev.off()
  pdf(file = paste0("quadrant_",group,".pdf"),paper = "USr")
  print(quadrant.list[[group]])
  dev.off()
}
ecdf.data <- na.omit(ecdf.data)
ecdf.data <- ecdf.data[is.finite(ecdf.data$data),]
ecdf.list[["combined"]] <- ECDF_plot(ecdf.data,ecdf.data$data,ecdf.data$group)
pdf(file = paste0("ecdf_","combined.pdf"),paper = "USr")
print(ecdf.list[["combined"]])
dev.off()
save(design.matrix,compare.list,heatmap_dm.list,heatmap_de.list,volcano_dm.list,ecdf.list,quadrant.list,file = output.Rdata)
