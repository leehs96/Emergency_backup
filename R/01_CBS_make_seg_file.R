library("DNAcopy")
args<-commandArgs(trailingOnly=TRUE)
#args[1]=input
#args[2]=sample_id
options("scipen"=100, "digits"=4)

print("start CBS")
cn <- read.table(args[1], header=F)
CNA.object <-CNA( genomdat = cn[,7], chrom = cn[,1],maploc = cn[,2], data.type = 'logratio', sampleid = args[2])
CNA.smoothed <- smooth.CNA(CNA.object)
segs <- segment(CNA.smoothed, verbose=2, min.width=2)
segs2 = segs$output
colnames(segs2) <- c("Sample", "Chromosome", "Start Position", "End Position","Num markers", "Seg.CN")
filename <- paste0(args[2],'.cnv.seg.tsv')
write.table(segs2[,1:6], file=filename, row.names=F, col.names=T, quote=F, sep="\t")

markername <- paste0(args[2],'.cnv.marker.tsv')
markers <- data.frame(paste(segs$data$chrom, segs$data$maploc, sep = ":"), segs$data$chrom, segs$data$maploc)
colnames(markers) <- c("Marker Name", "Chromosome", "Marker Position")
write.table(markers, file = markername, quote = FALSE, row.names = FALSE,sep = "\t")

print("Done!")
