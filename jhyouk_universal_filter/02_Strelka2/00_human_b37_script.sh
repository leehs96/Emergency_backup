normal=$1
tumor=$2
sampleName=$3
set -e

log=$3.strelka2.log

echo "strelka2 start for" $3 > $log
#/home/users/jhyouk/tools/Strelka2/strelka-2.9.9.centos6_x86_64/bin/configureStrelkaSomaticWorkflow.py --normalBam $1 --tumorBam $2 --ref /home/users/jhyouk/99_reference/human/GRCh37/human_g1k_v37.fasta --runDir $3 &>> $log

#for PCAWG reference genome including NC_007605 and hs37d5
/home/users/jhyouk/tools/Strelka2/strelka-2.9.9.centos6_x86_64/bin/configureStrelkaSomaticWorkflow.py --normalBam $1 --tumorBam $2 --ref /home/users/jhyouk/99_reference/human/hg19_PCAWG/genome.fa --runDir $3 &>> $log
echo "1st step done" >> $log
echo "runWorkflow" >> $log
$3/runWorkflow.py -m local -j 4 &>> $log
echo "strelka2 finish for" $3 >> $log
mv $3.strelka2.log $3.strelka2.success
