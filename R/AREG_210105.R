library(ggplot2)
library(reshape)
library(ggpubr)
library(rstatix)
library(dplyr)
library(RColorBrewer)
library(gplots)
library(pheatmap)

setwd('~/SNUH_data/')
areg <- read.csv('AREG_box.csv', header = T, row.names = 1)

areg <- as.data.frame(t(areg))
areg$AREG <- as.numeric(areg$AREG)

normal_box <- ggboxplot(
  areg,
  x = "Subtype",
  y = "AREG",
  fill = "Subtype") + 
  theme(
    legend.position = "none",
    plot.title = element_text(size = rel(1.5), hjust = 0.5),
    axis.title = element_text(size = rel(1.25))
  ) +
  ylab("AREG(FPKM)") +
  xlab("Subtype") +
  ggtitle("AREG expression")
normal_box + geom_jitter(alpha = 0.4, color = "#2e4057" , position=position_jitter(0.35))

stat.test_all_1 <- areg %>% t_test(AREG ~ Subtype, ref.group = "NT") 
stat.test_all_1 <- stat.test_all_1 %>% add_xy_position(x = "Subtype")
stat.test_all_1
normal_box + geom_jitter(alpha = 0.4, color = "#2e4057" , position=position_jitter(0.35)) + stat_pvalue_manual(stat.test_all_1, label = "p.adj.signif", tip.length = 0.01)



df <- read.csv('~/SNUH_data/merge_201228.csv', header = T)
#########################################################################
path <- "/users/hslee/Gene_set/jo/"
data1 <-list.files(path = path, pattern = '*.txt', full.names = T)
############################################
for (i in seq_along(data1)){
  a <- gsub(path,"",data1[i])
  b <- gsub('.txt','',a)
  name <- gsub('/','',b)
  assign(name, read.csv(data1[i], header = T, sep = ','))
}
############################################################



Heat <- function(df, gene_set, standard){
  ifelse(!(standard %in% gene_set[[1]]),gene_set[nrow(gene_set) + 1,] <- standard,NA)
  gene_df <- merge(df, gene_set, by= 1 )
  row.names(gene_df) <- gene_df$X
  gene_df <- gene_df %>% dplyr::select(!X)
  gene_df_t <- t(gene_df)
  gene_df_t <- as.data.frame(gene_df_t)
  gene_df_t_d <- gene_df_t %>% arrange(gene_df_t[names(gene_df_t)==standard])
  gene_df_t_d_t <- t(gene_df_t_d)
  gene_df_t_d_t_f <- gene_df_t_d_t[!rownames(gene_df_t_d_t)==standard,]
  gene_df_t_d_t_f <- filter(as.data.frame(gene_df_t_d_t_f), as.data.frame(gene_df_t_d_t_f) != 0)
  
  p <- pheatmap(
    gene_df_t_d_t_f,
    color = greenred(11),
    breaks = c(-1, -0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1),
    clustering_method = "ward.D2",
    border_color = NA,
    clustering_distance_rows = "correlation",
    fontsize_row = 10,
    cluster_rows = T,
    cluster_cols = F,
    scale = "row",
    show_colnames = F,
    show_rownames = T,
    legend = F,
  )
}


Heat(df, ADAM, flag)

cancer <- list(df)
names(cancer) <-c('SNUH_PTC,ADC')
flag <- c("AREG")
geneset <- list(ADAM,ALANINE_ASPARTATE_AND_GLUTAMATE_METABOLISM,angiogenesis,chemokine,efferocytosis,ERK,glycine_serine_threonine_metabolism,M2_gene_list,MAPK_pathway,MMP,PI3K_AKT,succci_goi)
names(geneset) <-c('ADAM','ALANINE_ASPARTATE_AND_GLUTAMATE_METABOLISM','angiogenesis','chemokine','efferocytosis','ERK','glycine_serine_threonine_metabolism','M2_gene_list','MAPK_pathway','MMP','PI3K_AKT','succci_goi')
############################################
setwd('~/조선욱_교수님/210105/')

for (i in 1:length(cancer)){
  for (j in 1:length(geneset)){
    png(filename=paste0(names(cancer)[[i]],'_',names(geneset)[[j]],'_',flag,'.png'),width=850,height=750,units = "px", bg="white" , res = 120,)
    Heat(cancer[[i]],geneset[[j]],flag)
    dev.off() 
  }
}
