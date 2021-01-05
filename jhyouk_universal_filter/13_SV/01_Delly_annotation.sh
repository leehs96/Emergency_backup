#!/bin/bash

# change reference fai file

sampleID=$1
tumor_col=$2
tumorBam=$3
normalBam=$4
pon=$5
scriptDir=$6
reference=$7 # mm10 or hg19

outDir=$(dirname $1)
rawDIR="~/cancer/WGS_hg37/variant_calling/delly_results/BCF"
log=$outDir/$1.SVprocessing.log

# START
echo $1 $3 $4 $5 > $log
#echo "start SV annotation from bcf combine"
#bcftools concat -a -O v -o $1.delly.vcf $rawDIR/$1.DEL.bcf $rawDIR/$1.DUP.bcf $rawDIR/$1.INS.bcf $rawDIR/$1.INV.bcf $rawDIR/$1.TRA.bcf &&
echo "done"
echo "filter somatic"
python $6/01.filter_somatic_delly.py $1.delly.vcf $2 && 
echo "done"
echo "Starting:sorting"
python $6/02.sorting_delly.py $1.delly.vcf.somatic &&
rm $1.delly.vcf.somatic
echo "done"
echo "Starting:annotate PON"
python $6/03.annotate_PON.py $1.delly.vcf.somatic.sort $5 ~/ref/hg37/human_g1k_v37.fasta.fai &&
rm $1.delly.vcf.somatic.sort
echo "done"
echo "Starting:find BP"
python $6/04.find_BP.py $1.delly.vcf.somatic.sort.pon $3 $4 &&
rm $1.delly.vcf.somatic.sort.pon
echo "done"
echo "Starting:BP adjustment"
python $6/05.BP_adjustment.py $1.delly.vcf.somatic.sort.pon.BPinfo &&
rm $1.delly.vcf.somatic.sort.pon.BPinfo
echo "done"
echo "Starting:count fragments and find newBP"
python $6/06.07.count_frag_find_newBP.py $1.delly.vcf.somatic.sort.pon.BPinfo.BPadj $3 $4 &&
rm $1.delly.vcf.somatic.sort.pon.BPinfo.BPadj
echo "done"
echo "Starting:annotate background information"
python $6/08.annotate_BGinfo.py $1.delly.vcf.somatic.sort.pon.BPinfo.BPadj.SVvaf $3 $4 &&
rm $1.delly.vcf.somatic.sort.pon.BPinfo.BPadj.SVvaf
echo "done"
echo "Starting:annotate paired normal same clipping"
python $6/09.PN_same_clipping.py $1.delly.vcf.somatic.sort.pon.BPinfo.BPadj.SVvaf.bginfo $3 $4 &&
rm $1.delly.vcf.somatic.sort.pon.BPinfo.BPadj.SVvaf.bginfo
echo "done"
echo "Starting:annotate mapq starting position variance"
python $6/10.MAPQ_start_pos.py $1.delly.vcf.somatic.sort.pon.BPinfo.BPadj.SVvaf.bginfo.pnsc $3 &&
rm $1.delly.vcf.somatic.sort.pon.BPinfo.BPadj.SVvaf.bginfo.pnsc
mv $1.delly.vcf.somatic.sort.pon.BPinfo.BPadj.SVvaf.bginfo.pnsc.mqpos $1.delly.vcf.somatic.annotated
echo "All Done"
####################여기까지함#########################################

echo "filter1 start"
python ~/script/jhyouk_universal_filter/13_SV/02_SV_filter1.py $1  &&
echo "filter1 done"
echo "06_invivo_column_Rearrange_afterfilter1"
python ~/script/jhyouk_universal_filter/13_SV/06_invivo_rearrange.py $1.delly.somatic.annotated.filter1.vcf &&
echo "06 done"


