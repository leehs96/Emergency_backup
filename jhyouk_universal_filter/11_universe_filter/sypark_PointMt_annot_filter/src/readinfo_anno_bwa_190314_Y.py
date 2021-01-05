
# coding: utf-8

# In[ ]:


#Arg1:point mutation vcf
#Arg2:tumor bam file

#for BWA
#updates
#180503: cigar loop update: error correction for deletion spanning the alt
#180518: cigar loop multiple error correction;add var_loca_reverse;change output format; get rlength from read
#180518: add informations: presence of clipping(soft of hard), INS, DEL in ref of var reads 
#180531: add informations: min, max location, refNM, varNM
#180717: indel is not available now ( write NA at all columns)
#180830: remove indel nt count from nm count
#180903: support gzip input
#181010: print NA when variant read count was zero.
#181022: indel can be annotated
#181101: unmapped read, duplicate read -> 'blank'
#181114: correction for multiple nucleotides in mutect2 such as CTT > CT, CT>AG
#181128: deletion length calculation error correction
#181129: deletion and insertion variant read detection error correction
#190221: count by readname

#by Youk
#190225: output_file name change (delete '.vcf')
#fetch(pos1-1,pos1+1) add for softclipping

import sys,pysam,gzip
from numpy import median
print(sys.argv[1])
if sys.argv[1][-3:]=='.gz':
    in_file=gzip.open(sys.argv[1])
else:
    in_file=open(sys.argv[1])
bam_file=pysam.AlignmentFile(sys.argv[2],'rb') #bam file
if sys.argv[1][-3:]=='.gz':
    out_file=open(sys.argv[1][:-3]+'.readinfo','w')
else:
    out_file=open(sys.argv[1]+'.readinfo','w')



in_line=in_file.readline().strip()
while in_line:
    if in_line[0:4]=='#CHR':
        out_file.write(in_line+'\ttemp_ref_readN\ttemp_var_readN\tref_medMQ\tvar_medMQ\tvar_LocaLeftMin;Med;Max\tvar_LocaRightMin;Med;Max\trefClip;pct;varClip;pct\trefIns;pct;varIns;pct\trefDel;pct;varDel;pct\tref_medMismatch\tvar_medMismatch\n')
    elif in_line[0]=='#':
        out_file.write(in_line+'\n')
    else:
        in_indi=in_line.split('\t')
        chr1=in_indi[0]
        pos1=int(in_indi[1])
        ref_nt=in_indi[3]
        alt_nt=(in_indi[4].split(',')[0]).split('/')[0]
        if len(ref_nt) == 1 and len(alt_nt)==1:
            mttype='snv'
        elif len(ref_nt) == len(alt_nt) and len(ref_nt) >1:
            mttype='snv'
            ref_nt=ref_nt[0]
            alt_nt=alt_nt[0]
        elif len(ref_nt) > len(alt_nt):
            mttype='del'
            cor_len=len(alt_nt)-1
            alt_nt=alt_nt[0]
            ref_nt=ref_nt[0:len(ref_nt)-cor_len]
            del_len=len(ref_nt)-1
        elif len(ref_nt) < len(alt_nt):
            mttype='ins'
            cor_len=len(ref_nt)-1
            ref_nt=ref_nt[0]
            alt_nt=alt_nt[0:len(alt_nt)-cor_len]
            ins_len=len(alt_nt)-1
            ins_seq=alt_nt[1:]
        else:
            print('ERROR: unknown mutation type. exiting')
            print(in_line)
            sys.exit(1)
        var_mapq=[];ref_mapq=[];var_loca_lt=[];var_loca_rt=[];var_nm=[];ref_nm=[]
        ref_n_list=[];var_n_list=[]
        ref_i=0;ref_d=0;ref_c=0;var_i=0;var_d=0;var_c=0
        for read in bam_file.fetch(chr1,pos1-1,pos1+1):
            if read.is_unmapped == True or read.is_duplicate == True: continue
            var_read='off'
            est_dist=pos1-read.reference_start-1
            rlength=read.infer_query_length()
            if read.cigartuples==None:
                continue
                
            if read.reference_start == pos1 and read.query_alignment_start == 0:
                continue              

            cigar_list=read.cigartuples
            current_m=0;current_i=0;current_d=0;target_del_stat=0
            for cigar in cigar_list:   #loop for calculate real distance 
                if cigar[0]==0 and (current_m + current_d)  <=est_dist:
                    current_m=current_m+cigar[1]
                elif cigar[0]==1 and (current_m + current_d) <=est_dist:
                    current_i=current_i+cigar[1]
                elif cigar[0]==2 and (current_m + current_d) <=est_dist:
                    if current_m+current_d+cigar[1] > est_dist:
                        target_del_stat=1
                        break
                    else:
                        current_d=current_d+cigar[1]
                elif current_m + current_d > est_dist:
                    break
                else:
                    'blank'
            if target_del_stat==1:
                ref_n_list.append(read.query_name)
            rel_dist=est_dist+current_i-current_d  # start with 0

            cigar_i=0;cigar_d=0;cigar_s=0;cigar_h=0; current_m=0;cigar_md=0  # cigar_md: any deletion after matching
            for cigar in cigar_list:   #loop for check presence of clipping, insertion, deletion
                if cigar[0]==0:
                    current_m += cigar[1]
                elif cigar[0]==1:
                    if mttype == 'ins' and pos1-read.reference_start==current_m+cigar_md and cigar[1]==ins_len and read.query_alignment_sequence[rel_dist+1:rel_dist+1+ins_len]==ins_seq:
                        var_read='on'
                    else:
                        cigar_i=cigar_i+1
                elif cigar[0]==2:
                    if mttype == 'del' and pos1-read.reference_start==current_m+cigar_md and cigar[1] == del_len:
                        var_read='on'
                    else:
                        cigar_d=cigar_d+1
                    if current_m > 0:
                        cigar_md += cigar[1]
                elif cigar[0]==4:
                    cigar_s=cigar_s+1
                elif cigar[0]==5:
                    cigar_h=cigar_h+1
            if target_del_stat !=1:
                try:
                    if mttype== 'snv' and read.query_alignment_sequence[rel_dist]==alt_nt:   #var_read
                        var_read='on'
                except:
                    continue		

            if var_read == 'on':
                var_n_list.append(read.query_name)
                var_mapq.append(read.mapping_quality)
                if read.has_tag('NM')==True:
                    current_nm=read.get_cigar_stats()[0][10] - (read.get_cigar_stats()[0][1]+read.get_cigar_stats()[0][2])
                    var_nm.append(current_nm)
                if cigar_list[0][0]==0:
                    dist=rel_dist
                    var_loca_lt.append(dist)
                    var_loca_rt.append(rlength-1-dist)
                elif cigar_list[0][0]==4 or cigar_list[0][0]==5:
                    dist=rel_dist+cigar_list[0][1]
                    var_loca_lt.append(dist)
                    var_loca_rt.append(rlength-1-dist)
                else:
                    'blank'

                if cigar_i > 0:
                    var_i = var_i+1
                if cigar_d > 0:
                    var_d = var_d+1
                if cigar_s > 0 or cigar_h > 0:
                    var_c = var_c+1
            if var_read == 'off':  #ref_read
                ref_n_list.append(read.query_name)
                ref_mapq.append(read.mapping_quality)
                if read.has_tag('NM')==True: 
                    current_nm=read.get_cigar_stats()[0][10] - (read.get_cigar_stats()[0][1]+read.get_cigar_stats()[0][2])
                    ref_nm.append(current_nm)
                if cigar_i > 0:
                    ref_i = ref_i+1
                if cigar_d > 0:
                    ref_d = ref_d+1
                if cigar_s > 0 or cigar_h > 0:
                    ref_c = ref_c+1
            else:
                'blank'
                
                
                
        if len(set(ref_n_list) & set(var_n_list)) > 0:  # in special case with different NT among overlapped pair
            var_n_list=list(set(var_n_list) - set(ref_n_list))
            ref_n_list=list(set(ref_n_list) - set(var_n_list))
        ref_n=len(list(set(ref_n_list)))
        var_n=len(list(set(var_n_list)))
        if len(ref_mapq)==0:
            mr_mq='NA'
        else:
            mr_mq=str(median(ref_mapq))
        if len(var_mapq)==0:
            mv_mq='NA'
        else:
            mv_mq=str(median(var_mapq))
        if len(var_loca_lt)==0:
            vlocal='NA'; vlocar='NA'
        else:
            vlocalmin=str(min(var_loca_lt))
            vlocalmed=str(median(var_loca_lt))
            vlocalmax=str(max(var_loca_lt))
            vlocarmin=str(min(var_loca_rt))
            vlocarmed=str(median(var_loca_rt))
            vlocarmax=str(max(var_loca_rt))
            vlocal=';'.join([vlocalmin, vlocalmed, vlocalmax])
            vlocar=';'.join([vlocarmin, vlocarmed, vlocarmax])
        if ref_n==0:
            ref_clip='NA;NA'
            ref_ins='NA;NA'
            ref_del='NA;NA'
        else:
            ref_clip=str(ref_c)+';'+str(round(ref_c*100/float(ref_n),2))
            ref_ins=str(ref_i)+';'+str(round(ref_i*100/float(ref_n),2))
            ref_del=str(ref_d)+';'+str(round(ref_d*100/float(ref_n),2))

        if var_n==0:
            var_clip='NA;NA'
            var_ins='NA;NA'
            var_del='NA;NA'
        else:
            var_clip=str(var_c)+';'+str(round(var_c*100/float(var_n),2))
            var_ins=str(var_i)+';'+str(round(var_i*100/float(var_n),2))
            var_del=str(var_d)+';'+str(round(var_d*100/float(var_n),2))
        clip_info=ref_clip+';'+var_clip
        ins_info=ref_ins+';'+var_ins
        del_info=ref_del+';'+var_del

        if len(ref_nm)==0:
            mr_nm='NA'
        else:
            mr_nm=str(median(ref_nm))
        if len(var_nm)==0:
            mv_nm='NA'
        else:
            mv_nm=str(median(var_nm))

        info_list=[str(ref_n), str(var_n), mr_mq, mv_mq, vlocal, vlocar, clip_info, ins_info, del_info, mr_nm, mv_nm ]
        out_file.write(in_line+'\t'+'\t'.join(info_list)+'\n')
    in_line=in_file.readline().strip()


