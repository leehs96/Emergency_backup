
# We just want every gene name with loci - not transcripts, etc.
wget ftp://ftp.ensembl.org/pub/grch37/current/gtf/homo_sapiens/Homo_sapiens.GRCh37.87.gtf.gz # get the file
gunzip Homo_sapiens.GRCh37.87.gtf.gz
awk '$3 == "gene"'Homo_sapiens.GRCh37.87.gtf > gene_only.gtf
# current convert2bed requires the transcript_id from a GTF file, but it can be blank
awk '{ if ($0 ~ "transcript_id") print $0; else print $0" transcript_id \"\";"; }' gene_only.gtf > gene2.gtf
convert2bed -i gtf < gene2.gtf > gene.bed
# now we just need the coords and gene name
awk -F '[\t;]' '{print $1,$2,$3,$12}' gene.bed > geneCleaned.bed 
# clean the file up - there are spaces as field separators
awk -F '[ ]' '{print $1,$2,$3,$6}' geneCleaned.bed > geneMinimal.bed
# remove quotes from gene names
sed 's/"//g' geneMinimal.bed > final.bed
