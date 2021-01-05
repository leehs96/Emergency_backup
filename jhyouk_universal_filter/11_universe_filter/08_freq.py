import sys
input_fn = sys.argv[1]
#input_fn = 'S1-2Gy-1.varscan.snp.vcf'
input_file = file(input_fn)
output_fn = input_fn.replace("_anno_v2.txt","_anno_add_vaf_v2.txt")
output_file = file(output_fn,'w')

input_line = input_file.readline().strip()
while input_line:
    if input_line[0:3] == 'Chr':
        #info = input_line+'\tfreq'
        info = 'Chr	Start	End	Ref	Alt	Func.refGene	Gene.refGene	GeneDetail.refGene	ExonicFunc.refGene	AAChange.refGene	cytoBand	SIFT_score	SIFT_pred	Polyphen2_HDIV_score	Polyphen2_HDIV_pred	Polyphen2_HVAR_score	Polyphen2_HVAR_pred	LRT_score	LRT_pred	MutationTaster_score	MutationTaster_pred	MutationAssessor_score	MutationAssessor_pred	FATHMM_score	FATHMM_pred	PROVEAN_score	PROVEAN_pred	VEST3_score	CADD_raw	CADD_phred	DANN_score	fathmm.MKL_coding_score	fathmm.MKL_coding_pred	MetaSVM_score	MetaSVM_pred	MetaLR_score	MetaLR_pred	integrated_fitCons_score	integrated_confidence_value	GERP.._RS	phyloP7way_vertebrate	phyloP20way_mammalian	phastCons7way_vertebrate	phastCons20way_mammalian	SiPhy_29way_logOdds	ExAC_ALL	ExAC_AFR	ExAC_AMR	ExAC_EAS	ExAC_FIN	ExAC_NFE	ExAC_OTH	ExAC_SAS	avsnp147	gerp..gt2	gnomAD_genome_ALL	gnomAD_genome_AFR	gnomAD_genome_AMR	gnomAD_genome_ASJ	gnomAD_genome_EAS	gnomAD_genome_FIN	gnomAD_genome_NFE	gnomAD_genome_OTH	Kaviar_AF	Kaviar_AC	Kaviar_AN	dann	dbscSNV_ADA_SCORE	dbscSNV_RF_SCORE	cosmic70	dpsi_max_tissue	dpsi_zscore	Otherinfo1	Otherinfo2	Otherinfo3	Otherinfo4	Otherinfo5	Otherinfo6	Otherinfo7	Otherinfo8	Otherinfo9	Otherinfo10	Otherinfo11	Otherinfo12	Otherinfo13	Otherinfo14	tumor_sample\tvaf'

        output_file.write(info+'\n')
    else:
        input_split = input_line.split('\t')
        input_vaf = input_split[85].split(':')[5]
        #print(input_vaf)
        vaf= '\t%s\n' %(input_vaf)
        output_file.write(input_line+vaf)
    input_line = input_file.readline().strip()
print('The End')
