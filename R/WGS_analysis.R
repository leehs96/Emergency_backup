path <- "/users/hslee/WGS_CE/"
data2 <-list.files(path = path, pattern = 'All_merged_anno_add_vaf_Cancer_genes.csv', full.names = T)
############################################
for (i in seq_along(data2)){
  a <- gsub(path,"",data2[i])
  b <- gsub('/','',a)
  name <- (gsub("_anno_add_vaf_Cancer_genes.csv","",b))
  a <- read.csv(data2[i], header = T , sep = ',')
  assign(name, a)
}




Cancer_gene1 <- read.csv('~/ref/Cancer_gene_data_COSMIC_201204.csv', sep = ',', header = T)
Cancer_gene <- select(Cancer_gene1,'Gene.Symbol', 'Tier',	'Hallmark', 'Somatic', 'Tumour.Types.Somatic.','Role.in.Cancer')
names(Cancer_gene)[names(Cancer_gene) == "Gene.Symbol"] <- "Gene.refGene"


All_merged1 <- merge(x= All_merged, y= Cancer_gene, by="Gene.refGene", all.x = T)
setwd('~/cancer/WGS_hg37/')
All_merged2 <- select(All_merged1, Chr, Start, Ref, Alt, Gene.refGene, Func.refGene, Role.in.Cancer, ExonicFunc.refGene,AAChange.refGene, gnomAD_genome_EAS, MetaLR_pred,dann, tumor_sample, t_vaf, Tier, Hallmark, Tumour.Types.Somatic.)

all_exonic <- subset(All_merged, All_merged$Func.refGene == 'exonic')
all_intronic <- subset(All_merged, All_merged$Func.refGene == 'intronic')
all_intergenic <- subset(All_merged, All_merged$Func.refGene == 'intergenic')

tail(names(sort(table(all_intronic$Gene.refGene))), 10)
tail(names(sort(table(all_exonic$Gene.refGene))), 10)
tail(names(sort(table(all_intergenic$Gene.refGene))), 10)


all_exonic_Can_gene <- filter(all_exonic, !is.na(Role.in.Cancer))
all_intronic_Can_gene <- filter(all_intronic, !is.na(Role.in.Cancer))

all_intronic_snp <- filter(all_intronic, Alt != '-' & Ref != '-')

length(which(duplicated(all_intronic_snp$Start)))



intron_snp1 <- read.table('~/WGS_CE/intron.prediction.snp_regSNP_intron.txt', sep = '\t', header = T)
intron_snp <- select(intron_snp1, chrom ,pos, alt, ref, disease, prob, tpr, fpr, splicing_site, strand)
intron_snp$chrom <- gsub('chr','', intron_snp$chrom)


names(intron_snp)[names(intron_snp) == "chrom"] <- "Chr"
names(intron_snp)[names(intron_snp) == "pos"] <- "Start"
names(intron_snp)[names(intron_snp) == "alt"] <- "Alt"
names(intron_snp)[names(intron_snp) == "ref"] <- "Ref"
names(intron_snp)[names(intron_snp) == "disease"] <- "regSNP_intron_disease"
names(intron_snp)[names(intron_snp) == "prob"] <- "regSNP_intron_Prob"
names(intron_snp)[names(intron_snp) == "tpr"] <- "regSNP_intron_TPR"
names(intron_snp)[names(intron_snp) == "fpr"] <- "regSNP_intron_FPR"
names(intron_snp)[names(intron_snp) == "splicing_site"] <- "regSNP_intron_splicing_site"
names(intron_snp)[names(intron_snp) == "strand"] <- "regSNP_intron_strand"




intron_snp$Func.refGene <- 'intronic'

All_merged3 <- merge(x = All_merged, y = intron_snp, by = c('Chr','Start','Alt', 'Ref','Func.refGene'), all.x = TRUE)

setwd('~/WGS_CE/')
write.csv(All_merged3,'All_merged_anno_add_vaf_Cancer_genes_regSNPintron.csv')



library(maftools)


path <- "/users/hslee/cancer/WGS_hg37/variant_calling/CNV_varscan/seg_file/"
data2 <-list.files(path = path, pattern = '*seg', full.names = T)
############################################
for (i in seq_along(data2)){
  a <- gsub(path,"",data2[i])
  b <- gsub('/','',a)
  name <- (gsub("varscan.cnv.seg","_seg",b))
  a <- read.csv(data2[i], header = T , sep = '\t',header = F)
  assign(name, a)ZZ
}

setwd('/users/hslee/cancer/WGS_hg37/variant_calling/CNV_varscan/')
colnames(name) <- c('name')


Uni <- read.csv('~/WGS_CE/All_merged_anno_add_vaf_Cancer_genes_regSNPintron.csv', sep = ',', header = T)

Exonic <- filter(Uni, Uni$Func.refGene == 'exonic')
Exonic_1 <- dplyr::select(Exonic,Chr,Start,End,Ref,Alt,Func.refGene,Gene.refGene,ExonicFunc.refGene,AAChange.refGene,t_vaf,tumor_sample,  SIFT_pred, Polyphen2_HDIV_pred, Polyphen2_HVAR_pred, LRT_pred, MutationTaster_pred, MutationAssessor_pred, 	FATHMM_pred, PROVEAN_pred, fathmm.MKL_coding_pred, MetaSVM_pred, MetaLR_pred, integrated_fitCons_score,  avsnp147, gerp..gt2, dann,	Kaviar_AN, gnomAD_genome_EAS, t_vaf, Role.in.Cancer, Tumour.Types.Somatic. )
