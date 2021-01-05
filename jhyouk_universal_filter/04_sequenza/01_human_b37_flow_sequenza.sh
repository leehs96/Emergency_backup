#!/bin/bash
set -e

mpileup_pwd=$1
sampleName=$2
normalPileup=$3

log=$3.sequenza.log

#create seqz file from the two mpileups
echo "bam2seqz for" $3 > $log
sequenza-utils bam2seqz -gc gc50base.hg37.wig.gz -n $1/$normalPileup.mpileup -t $1/$sampleName.mpileup -p -o $sampleName.seqz &>> $log

#binning
echo "binning" >> $log
sequenza-utils seqz_binning --window 100 --seqz $sampleName.seqz -o $sampleName.comp.seqz >> $log

# remove unwanted contigs in reference fasta
echo "mtglremove" >> $log
cat $sampleName.comp.seqz | grep -v 'MT' | grep -v 'GL' | grep -v 'JH' > $sampleName.comp.seqz.rmGLMTJH 2>> $log
#cat $sampleName.comp.seqz | grep -v 'chrM' | grep -v 'GL' | grep -v 'chrUn'| grep -v 'HLA' | grep -v 'KI' | grep -v 'JH' > $sampleName.comp.seqz.rmGLMTJH 2> $sampleName.mtglremove.out

## gzip 
echo "gzip" >> $log
gzip $sampleName.comp.seqz.rmGLMTJH >> $log

echo "done" >> $log

#cleanup
echo "cleanup" >> $log
rm $sampleName.seqz $sampleName.comp.seqz 

# further analysis
echo "Rscript" >> $log
#Rscript /home/users/jhyouk/06_mm10_SNUH_radiation/07_sequenza/03_Rscript_sequenza.R $sampleName.comp.seqz.rmGLMTJH.gz $sampleName &> $sampleName.Rscript.out
#Rscript /home/users/jhyouk/09_uveal_melanoma/04_sequenza/11_Rscript_sequenzawithpurity.R $sampleName.comp.seqz.rmGLMTJH.gz $sampleName-withpurity 1 2 &> $sampleName.purity.Rscript.out
Rscript /home/users/jhyouk/09_uveal_melanoma/04_sequenza/03_Rscript_sequenza.R $sampleName.comp.seqz.rmGLMTJH.gz $sampleName &>> $log
echo "Finish all for" $3 >> $log
mv $log $3.sequenza.success
