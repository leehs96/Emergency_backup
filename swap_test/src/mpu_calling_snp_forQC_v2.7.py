#python2.7
#last modification 2017-07-22 by SYP
import sys,os
import collections
vaf=0
depth=1
in_file=open(sys.argv[1])
out_file=open(sys.argv[1]+'.call','w')
out_file.write('#chr\tpos\tref\tvar\tdepth\tVAF\tref_forward\tref_reverse\tvar_forward\tvar_reverse\n')
in_line=in_file.readline().strip()
in_indi=in_line.split('\t')
sampleNo=1
def is_int(s):
	try:
		int(s)
		return True
	except:
		return False

while in_line:
	in_indi=in_line.split('\t')
	if in_indi[2]=='N':
		'blank'
	else:
		if int(in_indi[3])==0:
			out_file.write(in_indi[0]+'\t'+in_indi[1]+'\t'+in_indi[2]+'\t-\t0\tNA\t0\t0\t0\t0\n')
		if int(in_indi[3])>0:
			if in_indi[2]=='A':
				rnts='a'
			elif in_indi[2]=='C':
				rnts='c'
			elif in_indi[2]=='G':
				rnts='g'
			elif in_indi[2]=='T':
				rnts='t'
				
			if '^' in in_indi[4]:
				while '^' in in_indi[4]:
					start_idx=int(in_indi[4].index('^'))
					in_indi[4]=in_indi[4][:start_idx]+in_indi[4][start_idx+2:]

			if '+' in in_indi[4]:
				while '+' in in_indi[4]:
					ins_idx=int(in_indi[4].index('+'))
					try:
						n=0
						size_ins=0
						while is_int(in_indi[4][ins_idx+n+1])==True:
							n=n+1
						f=n-1
						for c in range(1,n+1):
							size_ins=size_ins+(int(in_indi[4][ins_idx+c])*(10**f))
							f=f-1	
						in_indi[4]=in_indi[4][:ins_idx]+in_indi[4][ins_idx+n+1+size_ins:]
						
					except Exception as e:
						print(e)
						sys.exit()		
	
			if '-' in in_indi[4]:
				while '-' in in_indi[4]:
					del_idx=int(in_indi[4].index('-'))
					try:
						n=0
						size_del=0
						while is_int(in_indi[4][del_idx+n+1])==True:
							n=n+1
						f=n-1
						for c in range(1,n+1):
							size_del=size_del+(int(in_indi[4][del_idx+c])*(10**f))
							f=f-1
						in_indi[4]=in_indi[4][:del_idx]+in_indi[4][del_idx+n+1+size_del:]
					except Exception as e:
						print(e)
						sys.exit()
	
			indi_nt=list(in_indi[4])
			num_nt=collections.Counter(indi_nt)
			ref_f=num_nt['.'];ref_r=num_nt[','];a_f=num_nt['A'];a_r=num_nt['a'];g_f=num_nt['G'];g_r=num_nt['g']
			t_f=num_nt['T'];t_r=num_nt['t'];c_f=num_nt['C'];c_r=num_nt['c']
			a_c=a_f+a_r;g_c=g_f+g_r;c_c=c_f+c_r;t_c=t_f+t_r
			dp=int(in_indi[3])
			lfc=[a_f,g_f,c_f,t_f]
			lrc=[a_r,g_r,c_r,t_r]
			lc=[a_c,g_c,c_c,t_c]
			ln=['A','G','C','T']
			if rnts=='a':
				ln[0]='.'
			elif rnts=='g':
				ln[1]='.'
			elif rnts=='c':
				ln[2]='.'
			elif rnts=='t':
				ln[3]='.'
			
			if int(dp)>=depth:  # Setting cut-off value of depth
				if a_c+g_c+c_c+t_c==0:
					out_file.write(in_indi[0]+'\t'+in_indi[1]+'\t'+in_indi[2]+'\t'+'-\t'+str(dp)+'\t'+str(round(0/float(dp),4))+'\t'+str(ref_f)+'\t'+str(ref_r)+'\t0\t0\n')
				else:
					for b in range (0,4):
						if ln[b]=='.' or lc[b]==0:
							'blank'
						else:
							if lc[b]/float(dp)>=vaf:   # Setting cut_off value of VAF
								out_file.write(in_indi[0]+'\t'+in_indi[1]+'\t'+in_indi[2]+'\t'+ln[b]+'\t'+str(dp)+'\t'+str(round(lc[b]/float(dp),4))+'\t'+str(ref_f)+'\t'+str(ref_r)+'\t'+str(lfc[b])+'\t'+str(lrc[b])+'\n')
							else:
								'blank'
	in_line=in_file.readline().strip()
