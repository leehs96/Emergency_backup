bam=$1
sampleID=$2

log=$sampleID.mpileup.log

echo "start mpileup for" $sampleID > $log

samtools mpileup -B -Q 20 -q 20 -f /home/users/jhyouk/99_reference/human/GRCh37/human_g1k_v37.fasta $1 -o $2.mpileup &>> $log
echo "finish mpileup for" $sampleID >> $log
mv $log $2.mpileup.success
