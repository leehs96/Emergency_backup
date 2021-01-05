
# coding: utf-8

# In[2]:

import sys
input_fn = sys.argv[1]
#input_fn = 'S1-2Gy-1.varscan.snp.vcf'
input_file = file(input_fn)
output_fn = input_fn.replace(".vcf",".filter2.vcf")
output_file = file(output_fn,'w')

input_line = input_file.readline()
while input_line:
    if input_line[0:2] == '##':
        'blank'
    elif input_line[0:6] == '#CHROM':
        output_file.write(input_line)
    else:
        input_split = input_line.split('\t')
        input_true = input_split[34]
        if input_true[0] == 'T':
            output_file.write(input_line)
        else:
            'blank'
    input_line = input_file.readline()
print('The End')
