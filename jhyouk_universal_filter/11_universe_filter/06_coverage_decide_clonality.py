
# coding: utf-8

# In[3]:


#Arg1 = vcf_filename
#Arg2 = tumor_bam
#Arg3 = mpileup.100kb.mpileup.100kbcov.covstat

def both_coverage(input_line,t_bam):
    input_split = input_line.split('\t')
    input_chr = input_split[0]
    input_pos = long(input_split[1])
    
    lt_coverage = []
    for pileupcolumn in t_bam.pileup(input_chr,input_pos-1000,input_pos,min_mapping_quality=0,min_base_quality=10):
        if input_pos-1000 <= pileupcolumn.pos and pileupcolumn.pos <= input_pos-1:
            lt_coverage.append(pileupcolumn.n)
            
    rt_coverage = []    
    for pileupcolumn in t_bam.pileup(input_chr,input_pos,input_pos+1000,min_mapping_quality=0,min_base_quality=10):
        if input_pos <= pileupcolumn.pos and pileupcolumn.pos <= input_pos+1000:
            rt_coverage.append(pileupcolumn.n)            

    return [numpy.mean(lt_coverage),numpy.mean(rt_coverage)]


import sys, pysam, gzip, numpy

input_fn = sys.argv[1]
#input_fn = '/home/users/jhyouk/06_mm10_SNUH_radiation/31_2_SNP_updated_190315/Liver_20Gy1_SO1_snp_union_2.readinfo.readc.rasmy_PanelofNormal.filter1.vcf'
input_file = file(input_fn)

t_bam_fn = sys.argv[2]
#t_bam_fn='/home/users/team_projects/Radiation_signature/02_bam/Liver_20Gy1_SO1.s.md.ir.br.bam'
t_bam = pysam.AlignmentFile(t_bam_fn,'rb')


info_fn = sys.argv[3]
#info_fn = "Liver_2Gy2_SO1.mpileup.100kbcov.covstat"
info_file = file("/home/users/jhyouk/06_mm10_SNUH_radiation/07_sequenza/"+info_fn)
info_line = info_file.readline().strip()
info_line = info_file.readline().strip()
info_coverage = float(info_line.split('\t')[3])

output_file = file(input_fn.replace('.vcf','')+'.coverage.vcf','w')

input_line = input_file.readline().strip()
prev_chr = '0'

j=0
while input_line[0:1] == '#':
    output_file.write(input_line + '\tT_Normalized_CN\tT_var_cell_portion\n')
    input_line = input_file.readline().strip()

clonal_cutoff=0.30 * 2 #Clonality (not vaf, corrected by copy number)

while input_line:
    input_split = input_line.split('\t')
    
    if input_split[0] != prev_chr:
        print input_split[0]
        prev_chr = input_split[0]
    
    #if input_split[1] != '3487989':
        #input_line = input_file.readline().strip()
        #continue

    info = '\tNA\tNA\n'
    if input_split[34] == 'F':
        output_file.write(input_line + info)
        
    else:
        input_tvaf = float(input_split[33])
        
        mean_coverage = both_coverage(input_line,t_bam)
        lt_mean = mean_coverage[0]*2 / info_coverage
        rt_mean = mean_coverage[1]*2 / info_coverage
        mean_cnv = max(round((lt_mean + rt_mean)/2),1)
        tvar_cor = input_tvaf * mean_cnv
                       
        
        info = '\t%s\t%s\n' % (mean_cnv, tvar_cor)
        output_file.write(input_line + info)
    input_line = input_file.readline().strip()

print 'THE END'

            

