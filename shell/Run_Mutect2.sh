tumor=$1
normal=$2
BAM_pass=$3

gatk Mutect2 -R ~/ref/hg37/human_g1k_v37.fasta -I $3/$1.s.md.br.bam -tumor $1 -I $3/$2.s.md.br.bam -normal $2  --germline-resource ~/ref/hg37/af-only-gnomad.raw.sites.vcf --panel-of-normals ~/ref/hg37/Mutect2-WGS-panel-b37.vcf -O ./$1.Mutect2.vcf
