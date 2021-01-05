tumor_bam=$1
normal_bam=$2
sampleID=$3

log = $sampleID.delly.log
echo "START DELLY for" $sampleID > $log

/home/users/tools/delly-0.7.6/delly/src/delly call -t DEL -q 15 -o $3.DEL.bcf -g /home/users/jhyouk/99_reference/human/hg19_PCAWG/genome.fa $1 $2 &> $3.DEL.out & 
/home/users/tools/delly-0.7.6/delly/src/delly call -t INS -q 15 -o $3.INS.bcf -g /home/users/jhyouk/99_reference/human/hg19_PCAWG/genome.fa $1 $2 &> $3.INS.out &
/home/users/tools/delly-0.7.6/delly/src/delly call -t DUP -q 15 -o $3.DUP.bcf -g /home/users/jhyouk/99_reference/human/hg19_PCAWG/genome.fa $1 $2 &> $3.DUP.out &
/home/users/tools/delly-0.7.6/delly/src/delly call -t INV -q 15 -o $3.INV.bcf -g /home/users/jhyouk/99_reference/human/hg19_PCAWG/genome.fa $1 $2 &> $3.INV.out &
/home/users/tools/delly-0.7.6/delly/src/delly call -t TRA -q 15 -o $3.TRA.bcf -g /home/users/jhyouk/99_reference/human/hg19_PCAWG/genome.fa $1 $2 &> $3.TRA.out
wait
echo 'Delly finish' $sampleID >> $log
mv $sampleID.delly.log $sampleID.delly.success
