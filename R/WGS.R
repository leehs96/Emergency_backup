library(ComplexHeatmap)
library(maftools)
library(dplyr)
library(reshape)
library(tximport)
library(RColorBrewer)
setwd('~/cancer/WGS_hg37/variant_calling/PON_filtered_vcf/InDel/for_anno/')
ADB_LL_indel <- read.table("ADB_LL.hg19_multianno.txt", sep = '\t', header = T)
var.annovar = system.file("extdata", "variants.hg19_multianno.txt", package = "maftools")
var.annovar.maf = annovarToMaf(annovar = var.annovar, Center = 'CSI-NUS', refBuild = 'hg19', 
                               tsbCol = 'Tumor_Sample_Barcode', table = 'ensGene')


plotmafSummary(maf = ADB_LL, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
getSampleSummary(var.annovar.maf)


t <- annovarToMaf(annovar = ADB_LL_indel, refBuild = 'hg19')


#path to TCGA LAML MAF file
laml.maf = system.file('extdata', 'tcga_laml.maf.gz' , ) 
#clinical information containing survival information and histology. This is optional
laml.clin = system.file('extdata', 'tcga_laml_annot.tsv', package = 'maftools') 

laml = read.maf(maf = laml.maf, clinicalData = laml.clin)
laml


t <- read.maf


###################################################################################################################



#ADB_LL.hg19_multianno.txt

path <- "/users/hslee/cancer/WGS_hg37/anno_result/"
data <-list.files(path = path, pattern = 'All_merged_anno_v2.txt', full.names = T)
############################################
for (i in seq_along(data)){
  a <- gsub(path,"",data[i])
  b <- gsub('/','',a)
  name <- gsub("_anno_v2.txt","",b)
  a <- annovarToMaf(annovar = data[i], refBuild = 'hg19', tsbCol = 'tumor_sample')
  a <- read.maf(a)
  assign(name, a)
}

######################################################################################################

path <- "/users/hslee/cancer/WGS_hg37/variant_calling/PON_filtered_vcf/SNP/for_anno/"
data2 <-list()
data2 <-list.files(path = path, pattern = '*.txt', full.names = T)
############################################
for (i in seq_along(data2)){
  a <- gsub(path,"",data2[i])
  b <- gsub('/','',a)
  name <- (gsub(".hg19_multianno.txt",".snp",b))
  a <- read.table(data2[i], header = T , sep = '\t')
  assign(name, a)
}

ADB_LL <- rbind(ADB_LL.indel, ADB_LL.snp)
ADB_LM <- rbind(ADB_LM.indel, ADB_LM.snp)
ADB_UL <- rbind(ADB_UL.indel, ADB_UL.snp)
ADB_UM <- rbind(ADB_UM.indel, ADB_UM.snp)
LMS_Cancer <- rbind(LMS_Cancer.indel, LMS_Cancer.snp)
LMS_Recur <- rbind(LMS_Recur.indel, LMS_Recur.snp)
LMS_LN_3R1 <- rbind(LMS_LN_3R1.indel, LMS_LN_3R1.snp)
LMS_LN_3L1 <- rbind(LMS_LN_3L1.indel, LMS_LN_3L1.snp)
SJH_Cen1 <- rbind(SJH_Cen1.indel, SJH_Cen1.snp)
SJH_Isth <- rbind(SJH_Isth.indel, SJH_Isth.snp)
SJH_Lat <- rbind(SJH_Lat.indel, SJH_Lat.snp)
SJH_LN <- rbind(SJH_LN.indel, SJH_LN.snp)




ADB <- rbind(ADB_LL,ADB_LM,ADB_UL,ADB_UM)


write.table(ADB, 'ADB_merged_anno_v2.txt', quote = F, row.names = F, sep = '\t')


SJH_LN <- SJH_LN[order(SJH_LN[,1], SJH_LN[, 2]),]





ADB_LL$tumor_sample <- 'ADB_LL'
ADB_LM$tumor_sample <- 'ADB_LM'
ADB_UL$tumor_sample <- 'ADB_UL'
ADB_UM$tumor_sample <- 'ADB_UM'
LMS_Cancer$tumor_sample <- 'LMS_Cancer'
LMS_Recur$tumor_sample <- 'LMS_Recur'
LMS_LN_3R1$tumor_sample <- 'LMS_LN_3R1'
LMS_LN_3L1$tumor_sample <- 'LMS_LN_3L1'
SJH_Cen1$tumor_sample <- 'SJH_Cen1'
SJH_Isth$tumor_sample <- 'SJH_Isth'
SJH_Lat$tumor_sample <- 'SJH_Lat'
SJH_LN$tumor_sample <- 'SJH_LN'

WGS_annovar <- rbind(ADB_LL,ADB_LM,ADB_UL,ADB_UM,LMS_Cancer,LMS_Recur,LMS_LN_3R1,LMS_LN_3L1,SJH_Cen1,SJH_Isth,SJH_Lat,SJH_LN)




dev.off()
setwd('~/cancer/WGS_hg37/')
write.table(WGS_annovar, 'All_merged_anno_v2.txt', quote = F, row.names = F, sep = '\t')

unique(WGS_annovar$tumor_sample)


plotmafSummary(maf = All_merged, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
laml.titv = titv(maf = All_merged, plot = T, useSyn = TRUE)


plotTiTv(res = laml.titv, showBarcodes = T)
#oncoplot for top ten mutated genes.
oncoplot(maf = All_merged, top = 20, SampleNamefontSize = 0.7, showTumorSampleBarcodes = T ,sampleOrder = c('ADB_LL','ADB_LM','ADB_UL','ADB_UM','LMS_Cancer','LMS_LN_3R1','LMS_LN_3L1','LMS_Recur','SJH_Cen1','SJH_Isth','SJH_Lat','SJH_LN'))
lollipopPlot(maf = All_merged, gene = 'AQP4', showDomainLabel = FALSE, AACol = "aaChange", labelPos = 104 )
laml.mutload = tcgaCompare(maf = All_merged, cohortName = 'CHA', logscale = TRUE, capture_size = 50)

rainfallPlot(maf = All_merged, detectChangePoints = TRUE, pointSize = 0.4, tsb = 'SJH_LN')

plotVaf(maf = All_merged, vafCol = 't_vaf')


#exclusive/co-occurance event analysis on top 10 mutated genes. 
somaticInteractions(maf = All_merged, top = 20, pvalue = c(0.05, 0.1))
laml.sig = oncodrive(maf = All_merged, AACol = 'aaChange', minMut = 1, pvalMethod = 'zscore')

plotOncodrive(res = laml.sig, fdrCutOff = 0.1, useFraction = F, labelSize = 0.2)
laml.pfam = pfamDomains(maf = All_merged, AACol = 'AAChange.refGene', top = 10)
laml.pfam$proteinSummary[,1:7, with = FALSE]

dgi = drugInteractions(maf = All_merged, fontSize = 0.75)
OncogenicPathways(maf = All_merged)
PlotOncogenicPathways(maf = All_merged, pathways = "RTK-RAS")
library("mclust")
ADB_LL = inferHeterogeneity(maf = All_merged, tsb = 'LMS_LN_3L1', vafCol = 'Otherinfo13')

library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg19, quietly = TRUE)


all.tnm = trinucleotideMatrix(maf = All_merged, ref_genome = 'BSgenome.Hsapiens.UCSC.hg19')

library(mclust)

ADB_LL.het = inferHeterogeneity(maf = All_merged, tsb = 'SJH_LN', vafCol = 't_vaf')
plotClusters(clusters = ADB_LL.het)


unique(WGS_annovar$cosmic70)









cancer/WGS_hg37/All_merged_anno_add_vaf_v2.txt

setwd('~/cancer/WGS_hg37/')
test <- read.table('All_merged_anno_add_vaf_v2.txt', sep = '\t', header = T)


#test$Chr <- paste0('chr', test$Chr)
test$vaf <- gsub('%', '', test$vaf)
test$vaf <- as.numeric(test$vaf)
test <- test %>% mutate(vafs = vaf / 100)

test <- select(test, -vaf)

#etwd()
write.table(test, "All_merged_anno_add_vaf_v2_chr.txt", sep = '\t', quote = F, row.names = F)

test$Chr <- paste0('chr', test$Chr) 

ADB_LL <- subset(test, test$tumor_sample == 'ADB_LL')
ADB_LM <- subset(test, test$tumor_sample == 'ADB_LM')
ADB_UL <- subset(test, test$tumor_sample == 'ADB_UL')
ADB_UM <- subset(test, test$tumor_sample == 'ADB_UM')
LMS_Cancer <- subset(test, test$tumor_sample == 'LMS_Cancer')
LMS_Recur <- subset(test, test$tumor_sample == 'LMS_Recur')
LMS_LN_3R1 <- subset(test, test$tumor_sample == 'LMS_LN_3R1')
LMS_LN_3L1 <- subset(test, test$tumor_sample == 'LMS_LN_3L1')
SJH_Cen1 <- subset(test, test$tumor_sample == 'SJH_Cen1')
SJH_Isth <- subset(test, test$tumor_sample == 'SJH_Isth')
SJH_Lat <- subset(test, test$tumor_sample == 'SJH_Lat')
SJH_LN <- subset(test, test$tumor_sample == 'SJH_LN')



write.table(SJH_LN, 'SJH_LN.anno.merged_chr.txt', sep = '\t', quote = F, row.names = F)


#test$freq <- as.numeric(test$freq)

all_arrange_vaf <- test %>% arrange(desc(freq))

all_exonic <- subset(test, test$Func.refGene == 'exonic')
all_intronic <- subset(test, test$Func.refGene == 'intronic')
all_intergenic <- subset(test, test$Func.refGene == 'intergenic')

library(ggpubr)


top10_intron <- list()
top10_intron <- tail(names(sort(table(all_intronic$Gene.refGene))), 30)

top10_intron_m <- subset(all_intronic, all_intronic$Gene.refGene%in%top10_intron)

(ggplot(top10_intron_m, aes(as.factor(Gene.refGene), fill=tumor_sample))
  + geom_bar()
  + geom_text(
    aes(label=''),
    stat='count'
  )  + scale_fill_brewer(palette = "Paired", direction = -1) + ggtitle("Top30 mutated gene \n Intronic") +
    xlab("Gene") + ylab("Counts") + theme(axis.text.x = element_text(angle = 45, hjust = 1))
) 


file_regSNP <- dplyr::select(all_intronic, Chr)


tail(names(sort(table(all_exonic$Gene.refGene))), 10)
tail(names(sort(table(all_intergenic$Gene.refGene))), 10)



intron_mut <- ggbarplot(all_intronic, x=)



BRAF <- subset(test, test$Gene.refGene == 'BRAF')

all_intronic$gnomAD_genome_EAS <- as.numeric(all_intronic$gnomAD_genome_EAS)
all_intronic_af_cutoff_0.001 <- dplyr::filter(all_intronic, all_intronic$gnomAD_genome_EAS < 0.0001 ) 

all_exonic_trans <- select(all_exonic, Chr, Start, Ref, Alt, Gene.refGene, ExonicFunc.refGene, gnomAD_genome_EAS, MetaLR_pred,dann, tumor_sample, t_vaf)




library(dplyr)



all_intronic_snv <-all_intronic %>% filter(all_intronic$Ref != '-' & all_intronic$Alt != '-') %>% select(Chr, Start, Ref, Alt, Gene.refGene,tumor_sample)
all_intronic_snv$Chr <- paste0('chr', all_intronic_snv$Chr)
setwd('~/cancer/WGS_hg37/')
write.csv(all_intronic_snv, 'intron_snv_file.csv')
length(which(duplicated(all_intronic_snv$Start)))



path <- "/users/hslee/cancer/WGS_hg37/variant_calling/CNV_varscan"
data2 <-list()
data2 <-list.files(path = path, pattern = '*ADB_LL.copynumber', full.names = T)
############################################
for (i in seq_along(data2)){
  a <- gsub(path,"",data2[i])
  b <- gsub('/','',a)
  name <- (gsub(".copynumber","_copy",b))
  a <- read.table(data2[i], header = T , sep = '\t')
  assign(name, a)
}




seg = system.file('extdata', 'SJH_LN.varscan.cnv.seg', package = 'maftools')
het = inferHeterogeneity(maf = All_merged, tsb = 'SJH_LN', segFile = seg, vafCol = 't_vaf')


#Visualizing results. Highlighting those variants on copynumber altered variants.
plotClusters(clusters = het, genes = 'CN_altered', showCNvars = T)
