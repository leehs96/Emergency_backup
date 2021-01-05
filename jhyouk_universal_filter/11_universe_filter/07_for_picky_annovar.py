import sys
input_fn = sys.argv[1]
#input_fn = 'S1-2Gy-1.varscan.snp.vcf'
input_file = file(input_fn)
output_fn = input_fn.replace("_union_2.readinfo.readc.rasmy_PanelofNormal.filter1.filter2.vcf",".for_annovar.vcf")
output_file = file(output_fn,'w')

input_line = input_file.readline()
while input_line:
    input_split = input_line.split('\t')
    CHROM = input_split[0]
    POS = input_split[1]
    ID = input_split[2]
    QUAL = input_split[3]
    REF = input_split[4]
    ALT = input_split[5]
    FILTER = 'PASS'
    INFO = '.'
    FORMAT = 'GT:GQ:DP:RD:AD:FREQ:DP4'
    info = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' %(CHROM,POS,ID,QUAL,REF,ALT,FILTER,INFO,FORMAT)
    output_file.write(info+'\n')
    input_line = input_file.readline()
print('The End')
