ID_1=$1
ID_2=$2
ref=$3
gtf=$4
outprefix=$5
STAR --genomeDir $3 --runThreadN 6 --sjdbGTFfile $4 --readFilesIn $1 $2  --sjdbOverhang 100 --outSAMtype BAM SortedByCoordinate --outFilterMultimapNmax 20 --outFilterType BySJout --twopassMode Basic  --outFileNamePrefix ./$5.STAR. > $5.star.out 
