import sys, itertools

depth_cut=10
cut_off1=int(sys.argv[2])

fn_file=open(sys.argv[1])  # id \t /path/to/input_file
print(sys.argv[1])
fn_line=fn_file.readline().strip()
sample_dics={}
sample_ids=[]
out_file=open('SampleSwapCheckResult.txt','w')
out_file.write('#Sample1\tSample2\tIdentity\tConsistency(cut-off='+str(cut_off1)+'%)\n')
print("Classifying SNPs")
while fn_line:
	fn_indi=fn_line.split('\t')
	sample_id=fn_indi[0]
	sample_ids.append(sample_id)
	list_ho=[];list_he=[];list_wt=[];list_total=[]
	in_file=open(fn_indi[1])
	in_line=in_file.readline().strip()
	while in_line:
		if in_line[0]=='#':
			'blank'
		else:
			in_indi=in_line.split('\t')
			idx='\t'.join(in_indi[0:4])
			if in_indi[5]=='NA':
				'blank'
			elif float(in_indi[4])<depth_cut:   # depth filtering
				'blank'
			else:
				list_total.append(idx)
				if in_indi[4]=='-':
					list_wt.append(idx)
				else:
					if float(in_indi[5])>=0.9:      # homozygote cut-off
						list_ho.append(idx)
					elif float(in_indi[5])>=0.1 and float(in_indi[5])<0.9: # heterozygote cut-off
						list_he.append(idx)
					elif float(in_indi[5])<0.1:     # wild-type cut-off
						list_wt.append(idx)
		in_line=in_file.readline().strip()
	sample_dics[sample_id]=[list_ho,list_he,list_wt,list_total]  ### order: ho, he, wt
	fn_line=fn_file.readline().strip()

sample_combi=list(itertools.combinations(sample_ids,2))
print("Calculating relationship")

for (a,b) in sample_combi:
	ho_ho=list(set(sample_dics[a][0]) & set(sample_dics[b][0]))
	fst_total=sample_dics[a][3]
	snd_total=sample_dics[b][3]
	identity=100*float(len(ho_ho))/(len(set(sample_dics[a][0]) & set(snd_total)))
	if identity>cut_off1:
		out_file.write(str(a)+'\t'+str(b)+'\t'+str(identity)+'%\tY\n')
	else:
		out_file.write(str(a)+'\t'+str(b)+'\t'+str(identity)+'%\tN\n')
		
		
