import subprocess,argparse,os
from multiprocessing import Pool


def run_mpu(ref, pos, sampleid, bamfile, srcdir):
	cmd='samtools mpileup -AB -q 20 -Q 20 -f '+ref+' -l '+pos+' -o '+sampleid+'.bqmpu '+bamfile+' && python '+srcdir+'/mpu_calling_snp_forQC_v2.7.py '+sampleid+'.bqmpu'
	p=subprocess.Popen(cmd, shell=True)
	p.wait()
	

def run_mpu_multicore(inputs, cores):
	pool = Pool(processes = cores)
	res = pool.starmap(run_mpu, inputs)
	


if __name__ == '__main__':
	parser=argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter, description=
		'this script check the sample swap using consistency of homozygote SNVs. Default loci of homozygote SNVs are derived from the result of 1000 genome project. Frequent variant loci were included.')
	parser.add_argument("BamList", type=argparse.FileType('r'), help='(Required) id\\tpath to input bam\\n should be written in the input text file line by line.(Example,\nID1\t/path/to/BAM1\nID2\t/path/to/BAM2\n...\t...\n)')
	parser.add_argument("-p","--PositionList", type=str, help="Position list for evaluation. Default('/users/hslee/script/swap_test/src/hg19_aeg.sites.2015_08_chrall_snp_cds_cut_0.1_0.9.txt.fi.segdup.pos')", default='/users/hslee/script/swap_test/src/hg19_aeg.sites.2015_08_chrall_snp_cds_cut_0.1_0.9.txt.fi.segdup.pos')
	parser.add_argument("-r","--reference", type=str, help="Reference fasta which was used for alignment. Default('/users/hslee/ref/hg37/human_g1k_v37.fasta')", default='/users/hslee/ref/hg37/human_g1k_v37.fasta')
	parser.add_argument("-s","--ScriptDir",type=str,help="The location of scripts for this tool. Default('/users/hslee/script/swap_test/src')", default='/users/hslee/script/swap_test/src')
	parser.add_argument("-c","--cutoff", type=int, help="Cutoff of percentage for detecting consistent samples. Default(90)", default=90)
	parser.add_argument("-t","--thread", type=int, help="Number of thread for parallel pileup. Default(1)", default=1)
	args=parser.parse_args()
	cwd=os.getcwd()
	n=0
	output_list=[]
	out_file=open('input_list_for_swap_check.tsv','w')
	input_list=[]
	in_line=args.BamList.readline().strip()
	while in_line:
		in_indi=in_line.split('\t')
		sampleid=in_indi[0]
		bamfile=in_indi[1]
		input_list.append([args.reference, args.PositionList, sampleid, bamfile, args.ScriptDir])
		out_file.write(sampleid+'\t'+sampleid+'.bqmpu.call\n')
		in_line=args.BamList.readline().strip()
	out_file.close()

	run_mpu_multicore(input_list, args.thread)
	p=subprocess.Popen('python '+args.ScriptDir+'/consistency_calc.py '+cwd+'/input_list_for_swap_check.tsv '+str(args.cutoff), shell =True)
	p.wait()
	
	
