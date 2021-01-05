#!/bin/bash
set -e

sampleID=$1
mttype=$2 # snp or indel
reference=$3 # mm10 or hg19
somaticbam=$4
germlinebam=$5
panelofnormal=$6

case $3 in
	mm10)	species="mouse";;
	hg19)	species="human";;
esac
case $2 in
	snp)	temp_mttype="snvs";;
	indel)	temp_mttype="indels";;
esac


outDir=$(dirname $1)
log=$outDir/$1.$2.annot.log


echo $1 $2 $3 $4 $5 > $log
echo varscan somatic filtering >> $log
python ~/script/jhyouk_universal_filter/11_universe_filter/00_varscan_somaticfilter.py ~/cancer/WGS_hg37/variant_calling/filtered/$1/$1.varscan.$2.pass.vcf &&
echo 'done' >> $log
echo "start: union of pass call in varscan2 somatic & strelka2" >>$log
python ~/script/jhyouk_universal_filter/11_universe_filter/00_vcf_combination_by_Youk_$2.py $1 $species 3 ~/cancer/WGS_hg37/variant_calling/filtered/$1/$1.strelka.somatic.$temp_mttype.pass.vcf ~/cancer/WGS_hg37/variant_calling/filtered/$1/$1.varscan.$2.pass.somatic.vcf ~/cancer/WGS_hg37/variant_calling/PON_filtered_vcf/Mutect2/pass/SNP/$1.Mutect2.pass.snp.vcf &&
echo "done" >>$log
echo "initial annotation" >> $log
sh ~/script/jhyouk_universal_filter/11_universe_filter/sypark_PointMt_annot_filter/PointMt_annot.sh $1_$2_union_3.vcf $4 $5 ~/script/jhyouk_universal_filter/11_universe_filter/sypark_PointMt_annot_filter/src $3 &&
echo "done" >>$log
echo "panel of normal annotation" >> $log
python2.7 ~/script/jhyouk_universal_filter/11_universe_filter/02_AddNpanelToVCF_$2.py $1_$2_union_3.readinfo.readc.rasmy.vcf $6 PanelofNormal $species &&
echo "done" >>$log
echo "filter1 using sample and germline information"  >>$log
python2.7 ~/script/jhyouk_universal_filter/11_universe_filter/03_$2_filter1.py $1_$2_union_3.readinfo.readc.rasmy_PanelofNormal.vcf &&
echo "filter1 done" >> $log
echo "remove intermediate files" >> $log
rm $1_$2_union_3.vcf
rm $1_$2_union_3.readinfo.readc.rasmy.vcf
rm $1_$2_union_3.readinfo.readc.rasmy_PanelofNormal.vcf
echo 'filter only PASS' >> $log
python2.7 ~/script/jhyouk_universal_filter/11_universe_filter/07_only_pass.py $1_$2_union_3.readinfo.readc.rasmy_PanelofNormal.filter1.vcf &&
echo 'make Annovar file' >> $log
python2.7 ~/script/jhyouk_universal_filter/11_universe_filter/07_for_picky_annovar.py $1_$2_union_3.readinfo.readc.rasmy_PanelofNormal.filter1.filter2.vcf &&
rm $1_$2_union_3.readinfo.readc.rasmy_PanelofNormal.filter1.vcf
echo "Finish!!" >> $log
