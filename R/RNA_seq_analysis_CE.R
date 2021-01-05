library(DESeq2)
library(ggrepel)
library(tximport)
library(RColorBrewer)
library(pheatmap)

setwd('/users/hslee/cancer/RNA_hg37/trim_qc/rsem/rsem_results')

files <- c(#"ADB_LL_.genes.results",
           #"ADB_LM_.genes.results",
           #"ADB_UL_.genes.results",
           #"ADB_UM_.genes.results",
           "ADB_NT_.genes.results",
           #"LMS_Cancer_.genes.results",
           #"LMS_LN_3R1_.genes.results",
           #"LMS_LN_3L1_.genes.results",
           #"LMS_Recur_.genes.results",
           "LMS_NT_.genes.results",
           "SJH_Cen1_.genes.results",
           "SJH_Lat_.genes.results",
           "SJH_Isth_.genes.results",
           "SJH_LN_.genes.results",
           "SJH_NT_.genes.results")

names(files) <- c(#"ADB_LL","ADB_LM","ADB_UL","ADB_UM",
                  "ADB_NT",
                  #"LMS_Cancer","LMS_LN_3R1","LMS_LN_3L1","LMS_Recur",
                  "LMS_NT",
                  "SJH_Cen1","SJH_Lat","SJH_Isth","SJH_LN",
                  "SJH_NT")

#tximport datafiles in directory
txData <- tximport(files,'rsem', txIn = T, txOut = T)

#change seqLen 0s to 1s
txData$length[txData$length == 0] <-1

#read column data
columns <- read.csv("column_SJH.csv", header = T, sep = "," , row.names = 1)

#data structure
dds <- DESeqDataSetFromTximport(txData, columns , ~Tissue)

#remove lowly exp. gene
keep <- rowSums(counts(dds)) >= 10

dds <- dds[keep,]

#export normalized exp values
ddsN <- estimateSizeFactors(dds)
ddsN <- estimateDispersions(ddsN)

Ndds <- counts(ddsN, normalized = T)

write.csv(Ndds,"CE_RNA_gene_norm_SJH.csv")

#DEG
ddsDE <- DESeq(dds)

#results change p-adj sig value 0.05
res <- results(ddsDE, contrast = c("Tissue","Cancer", "Normal") , alpha = 0.05)


#order by p-adj
resOrdered <- res[order(res$padj),]

write.csv(resOrdered,"desep_results_gene_SJH.csv")

##Plotting

#plot mean exp by lfc
plotMA(ddsDE)

#PCA
vsd <- vst(dds, blind = F)
plotPCA(vsd, intgroup = c("ID", "Tissue"))

? plotP
#sample to sample dist.
sampleDist <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDist)
rownames(sampleDistMatrix) <- vsd$Tissue
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")))(255)
#heat-map sample sample distance
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDist,
         clustering_distance_cols = sampleDist,
         cluster_cols = T,
         cluster_rows = T,
         border_color = "white",
         display_numbers = T,
         cutree_rows = 2,
         cutree_cols = 2,
         color = colors,
         show_rownames = T,
         show_colnames = T
)

#ggplot

library(ggplot2)
library(reshape)

setwd('~/cancer/RNA_hg37/trim_qc/rsem/rsem_results/')
##exp file
#normExp <- read.csv("CE_RNA_gene_norm_SJH.csv", row.names = 1)
normExp <- read.csv("CE_RNA_gene_norm_ADB.csv", row.names = 1)

##DESeq
#statsDE <- read.csv("desep_results_gene_SJH.csv",row.names = 1)
statsDE <- read.csv("desep_results_gene_SJH.csv",row.names = 1)
statsDE <- na.omit(statsDE)
statsDE <- statsDE[order(statsDE$padj),]
statsDE$sig <- ifelse(statsDE$padj <= 0.05, "sig", "not")


#plotMA replicate
plotMA <- ggplot(statsDE, aes(x = log10(baseMean) , y = log2FoldChange, color = sig)) + 
  geom_point() +
  theme(legend.position = "none")
plotMA

#Volcano
statsDE$genelabels <- rownames(statsDE) %in% rownames(statsDE[1:20,])

#set threshold
padj.cutoff <- 0.01
lfc.cutoff <- 0.58

threshold <- statsDE$padj < padj.cutoff & abs(statsDE$log2FoldChange) > lfc.cutoff & statsDE$baseMean >= 1
length(which(threshold))
statsDE$threshold <- threshold

volcano <- 
  ggplot(statsDE, aes(x = log2FoldChange, y = -log10(padj), color = threshold)) + 
  geom_point() +
  geom_text_repel(
    colour = 'black',
    size = 3,
    aes(
      x = log2FoldChange,
      y = -log10(padj),
      label = ifelse(genelabels == T, rownames(statsDE),""))
  ) +
  ggtitle("Cancer vs Normal") +
  xlab("Log2 fold change") +
  ylab("-Log10 adjusted p-value") +
  theme(
    legend.position = "none",
    plot.title = element_text(size = rel(1.5), hjust = 0.5),
    axis.title = element_text(size = rel(1.25))
  ) #+
  #scale_x_continuous(limits=c(-20,20), breaks=seq(-20,20,5))

volcano

