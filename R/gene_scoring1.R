library(RColorBrewer)
library(gplots)
library(ggplot2)
library(dplyr)
library(pheatmap)
library(ggpubr)
library(reshape)


df <- read.csv('~/TCGA_RNA_data/201228/etc/171013_FPKM.csv', header = T)
info <- read.csv('~/TCGA_RNA_data/201228/etc/info_201228.csv', header = T, row.names = 1)
df <- df[-which(duplicated(df$X)),]
rownames(df) <- df$X
df <- select(df, -X)


tt <- rbind(info, df1)

setwd('~/TCGA_RNA_data/201228/etc/')
write.csv(tt,'merge_201228.csv')
write.csv(info,'info.csv')

deg_gs <- read.csv('DEG.csv', header = T, row.names = 1)
gs <- read.csv('201228_gs.csv', header = T, row.names = 1)
annotation <- data.frame(sampletype = rep(c("Normal","DTC","ATC"),c(68,168,10)))
row.names(annotation) <- names(gs)
gs <- gs[rownames(gs)!='DEG',]

p <- pheatmap(
  gs,
  color = greenred(11),
  #annotation_col = annotation,
  breaks = c(-2,-0.8,-0.6,-0.3,-0.2,0,0.2,0.3,0.6,0.8,2),
  clustering_method = "ward.D2",
  border_color = NA,
  clustering_distance_rows = "correlation",
  fontsize_row = 12,
  cluster_rows = T,
  cluster_cols = F,
  scale = "row",
  show_colnames = F,
  show_rownames = T,
  legend = F,
  cutree_cols = 1,
  cutree_rows = 3
  )
p

