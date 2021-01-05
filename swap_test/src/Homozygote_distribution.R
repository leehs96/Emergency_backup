#!/home/users/kjyi/tools/R/R-3.4.0/bin/Rscript
require(tidyverse)
args <- commandArgs(trailingOnly = T)
fn_list <- read_tsv(args[1], col_names = F)  #mpileup.call list
id_list <- read_tsv(args[2], col_names = F)  #id list
colnames(fn_list) <- c("files")
colnames(id_list) <- c("ids")
p <- list()
for ( n in 1:nrow(fn_list)){
	print(n)
	dt <- read_tsv(fn_list$files[n])
	id <- id_list$ids[n]
	print(id)
	dt %>% filter(VAF > 0.5) -> subdt
	png(paste0(id,".homovaf.png"), width = 6, height = 4, units = "in", res = 600)
	g <- ggplot(subdt)+
		geom_histogram(aes(VAF*100), breaks = seq(50, 100, by = 0.1), closed = "right") +
		ggtitle(id) + 
		xlab("VAF (%)") +
		ylab("Number of SNPs")
	print(g)
	dev.off()
	}
