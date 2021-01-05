sampleID=$1
germlineID=$2
bam_pwd=$3
panelofnormal=$4
ref=$5

sh /home/users/jhyouk/81_filter_test_LADC/13_SV/01_Delly_annotation.sh $1 10 $bam_pwd/$1.s.md.ir.br.bam $bam_pwd/$2.s.md.ir.br.bam $panelofnormal /home/users/jhyouk/81_filter_test_LADC/13_SV/Delly_annotation_scripts $ref
