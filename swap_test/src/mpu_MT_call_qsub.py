import sys, os
in_file=open(sys.argv[1]) #bam file list (absolute path)
id_file=open(sys.argv[2]) #id list
in_line=in_file.readline().strip()
id_line=id_file.readline().strip()
reffa='/home/users/data/01_reference/human_g1k_v37/human_g1k_v37.fasta'
while in_line:
	out_file=open('mpu_MT_and_call_%s.sh' % id_line,'w')
	out_file.write('#!/bin/bash\n')
	out_file.write('#PBS -l nodes=1:ppn=1\n')
	out_file.write('#PBS -j oe\n')
	out_file.write('#PBS -o ./%s.MT.oe\n' % id_line)
	out_file.write('samtools mpileup -B -q 20 -Q 20 -r MT -f %s -o %s.MT.mpileup %s\n' % (reffa, id_line, in_line))
	out_file.write('python2.7 /home/users/sypark/01_Python_files/bamqc/mpu_calling_snp_forQC_v2.7.py %s.MT.mpileup$\n' % id_line)
	out_file.close()
	os.system('qsub -q week mpu_MT_and_call_%s.sh' % id_line )
	in_line=in_file.readline().strip()
	id_line=id_file.readline().strip()


