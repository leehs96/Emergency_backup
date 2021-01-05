
# coding: utf-8

# In[1]:


# input_fn : somatic + germline bam file name list (ex>input_germlineLiver_20Gy_1_G.s.md.ir.br.bam) delimited by tab
import sys
input_fn = sys.argv[1]
input_date = sys.argv[2] # 6digits ex>181120 for making run_script file
#input_fn = "181120_bam.txt"
#input_date = "181120"
input_file = file(input_fn)
#output_04 = file("04_strelka/01_2_run_strelka_" + input_date + "_human.sh",'w')
#output_05 = file("05_mutect/00_1_run_mutect_" + input_date + "_human.sh",'w')
output_03 = file("03_mpileup/00_1_run_mpileup_" + input_date + "_human.sh",'w')
output_04 = file("04_sequenza/01_1_run_sequenza_" + input_date + "_human.sh",'w')
output_01 = file("01_delly/01_run_Delly_" + input_date + "_human.sh",'w')
output_05 = file("05_varscan/00_1_run_varscan_" + input_date + "_human.sh",'w')
output_02 = file("02_Strelka2/00_1_run_strelka2_" + input_date + "_human.sh",'w')
#output_15 = file("15_MuTect2/00_1_run_Mutect2_"+ input_date + "_human.sh",'w')

input_line = input_file.readline().strip()
while input_line:
    input_somatic_bam = input_line.split('\t')[0]
    input_somatic = input_somatic_bam.split('/')[-1].split('.bam')[0]
    input_germline_bam = input_line.split('\t')[1]
    input_germline = input_germline_bam.split('/')[-1].split('.bam')[0]
    
    #output_04.write("sh 00_human_b37_strelka.sh /home/users/team_projects/uveal_melanoma/05_bam/" + input_germline_bam +                     " /home/users/team_projects/uveal_melanoma/05_bam/" + input_somatic_bam + " " + input_somatic +                     " /home/users/jhyouk/99_reference/human/hg38/Homo_sapiens_assembly38.fasta" + '\n')
    #output_05.write("sh 00_human_b37_Mutect.sh /home/users/team_projects/uveal_melanoma/05_bam " + input_somatic + " " + input_germline + '\n')
    output_03.write("sh 00_human_b37_script_mpileup.sh " + input_somatic_bam + ' ' +  input_somatic + '\n')  # Need to add new_germline
    output_03.write("sh 00_human_b37_script_mpileup.sh " + input_germline_bam + ' ' + input_germline + '\n')  # Need to add new_germline
    output_04.write("sh 01_human_b37_flow_sequenza.sh ../03_mpileup " + input_somatic + " " + input_germline + '\n')
    output_01.write("sh 00_human_b37_Delly.sh %s %s %s\n" %(input_somatic_bam, input_germline_bam, input_somatic))
    output_05.write("sh 00_human_b37_varscan.sh ../03_mpileup "+ input_somatic + " " + input_germline + '\n')
    output_02.write("sh 00_human_b37_script.sh %s %s %s\n" %(input_germline_bam, input_somatic_bam, input_somatic))
    #output_15.write("sh 00_human_b37_script_Mutect2.sh " + input_somatic + " " + input_germline + '\n')
    
    input_line = input_file.readline().strip()


print 'The End'
    
    

