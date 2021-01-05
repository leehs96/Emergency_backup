library(biomaRt)
library(dplyr)
setwd('~/cancer/RNA_hg37/trim_qc/rsem/rsem_results/')

df <- read.csv('desep_results_gene_ADB.csv')


human = useEnsembl(biomart = "ensembl", dataset="hsapiens_gene_ensembl") 
mouse = useEnsembl(biomart = "ensembl", dataset="mmusculus_gene_ensembl") 


gene <- getLDS(
   attributes = c("ensembl_gene_id"),
   filters = "ensembl_gene_id",
   values = df$X,
   mart = human,
   attributesL = c("hgnc_symbol"),
   martL = human,
   uniqueRows=T
)
#listAttributes(human)

#gene <- arrange(gene,Gene.stable.ID)


t <- merge(df, gene, by=1, all.x = T)

t1 <- na.omit(t)
t1 <- t1[,-1]
t1 <- t[order(t$padj),]
t1 <- t1 %>% relocate(HGNC.symbol)

write.csv(t1,'desep_results_gene_ADB(symbol).csv',row.names = F)



length(which(duplicated(t1$X)))




conversion <- function(df, species_c, species_a, current, alteration){
   
   ifelse(species_c=='human',data="hsapiens_gene_ensembl",data1="mmusculus_gene_ensembl")
   ifelse(species_a=='human',data="hsapiens_gene_ensembl",data2="mmusculus_gene_ensembl")
   
   gene <- getLDS(
      attributes = c(current),
      filters = current,
      values = df[,1],
      mart = useEnsembl(biomart = "ensembl", dataset=data1),
      attributesL = c(alteration),
      martL = useEnsembl(biomart = "ensembl", dataset=data2),
      uniqueRows=T
   )
   
   
}





###################################################################################
#ensembl <- useEnsembl(biomart = "ensembl", dataset="hsapiens_gene_ensembl") 
#chr1genes <- getBM(attributes = c('ensembl_gene_id','ensembl_transcript_id','hgnc_symbol',
#                                  'chromosome_name','start_position','end_position'), 
#                   filters ='chromosome_name', values ="1", mart = ensembl)
#
#######################################################################################
#humanx <- unique(genes[, 2])

