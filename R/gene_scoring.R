path <- "/users/hslee/TCGA_RNA_data/jo/"
path <- "/users/hslee/조선욱_교수님/201222/score/out/"
data2 <-list.files(path = path, pattern = '*TNBC.csv', full.names = T)
############################################
for (i in seq_along(data2)){
  a <- gsub(path,"",data2[i])
  b <- gsub('/','',a)
  name <- (gsub(".csv","_gs",b))
  a <- read.csv(data2[i], header = T, sep = ',', row.names = 1)
  a <- dplyr::select(a, grep('.01$', names(a)))
  #names(a) <- gsub('.01','',names(a))
  #a <- select(as.data.frame(t(a)), CXCL16, IL8, AREG)
  as.data.frame(assign(name, a))
}

for (i in seq_along(data2)){
  a <- gsub(path,"",data2[i])
  b <- gsub('/','',a)
  name <- (gsub("_RNAseq_TCGA.csv","",b))
  a <- read.csv(data2[i], header = T, sep = ',', row.names = 1)
  gene <- a[,1] 
  as.data.frame(assign(name, a))
}
