#!/bin/bash
input=$1
tbam=$2
nbam=$3
srcDir=$4
ref=$5

outDir=$(dirname $input)
log=$outDir/$input.annot.log

echo $input > $log
echo "Starting: Readinfo annot" >> $log
(python2.7 $srcDir/readinfo_anno_bwa_190314_Y.py $input $tbam) &&
echo "done" >> $log
echo "Starting: N readcount annot" >> $log
(python2.7 $srcDir/readcount_only_anno_190314_clipinsdeladd_Y.py $input.readinfo $nbam pairN) &&
echo "done" >> $log
#rm $input.readinfo
echo "Starting: Local reassembly annot" >> $log
(python2.7 $srcDir/read_local_reassembly_RPsize_190311.py $input.readinfo.readc $tbam $nbam $ref) &&
echo "done" >> $log
#rm $input.readinfo.readc

