#!/bin/bash
#PBS -l nodes=1:ppn=1
#PBS -j oe
#PBS -o /home/users/jhyouk/81_filter_test_LADC/10_PanalOfNormal/qsub_sdout/run.korean36.mpu_19.sh.sdout
cd /home/users/jhyouk/81_filter_test_LADC/10_PanalOfNormal
(samtools mpileup -AB -d 3000 -q 0 -Q 0 -f /home/users/jhyouk/99_reference/human/GRCh37/human_g1k_v37.fasta -r 20 /home/users/team_projects/BRCA_yonsei/02_bam/4B2.s.md.ir.br.bam /home/users/team_projects/BRCA_yonsei/02_bam/ON015.s.md.ir.br.bam /home/users/team_projects/BRCA_yonsei/02_bam/2B1.s.md.ir.br.bam /home/users/team_projects/uveal_melanoma/02_bam/1-17-Blood.s.md.ir.br.bam /home/users/team_projects/uveal_melanoma/02_bam/2-24-Blood.s.md.ir.br.bam /home/users/team_projects/uveal_melanoma/02_bam/3-29-Blood.s.md.ir.br.bam /home/users/team_projects/uveal_melanoma/02_bam/4-42-Blood.s.md.ir.br.bam /home/users/team_projects/uveal_melanoma/02_bam/5-68-Normal.s.md.ir.br.bam /home/users/team_projects/uveal_melanoma/02_bam/6-88-Blood.s.md.ir.br.bam /home/users/team_projects/uveal_melanoma/02_bam/8802875-Blood.s.md.ir.br.bam /home/users/team_projects/uveal_melanoma/02_bam/8769056-Blood.s.md.ir.br.bam /home/users/team_projects/uveal_melanoma/02_bam/8804915-Blood.s.md.ir.br.bam /home/users/team_projects/uveal_melanoma/02_bam/8808738-Blood.s.md.ir.br.bam /home/users/team_projects/uveal_melanoma/02_bam/8812875-Blood.s.md.ir.br.bam /home/users/team_projects/Radiation_signature/05_bam_human_sample_GRCh37/HBIR1_germline.s.md.ir.br.bam /home/users/team_projects/Radiation_signature/05_bam_human_sample_GRCh37/HBIR2_germ.s.md.ir.br.bam /home/users/team_projects/Radiation_signature/05_bam_human_sample_GRCh37/HCIR1.s.md.ir.br.bam /home/users/team_projects/Radiation_signature/05_bam_human_sample_GRCh37/hs_study5_germline_HCIR2.s.md.ir.br.bam /home/users/team_projects/Radiation_signature/05_bam_human_sample_GRCh37/hs_study5_colon_control_HCNor2_bulktissue.s.md.ir.br.bam /home/users/team_projects/uveal_melanoma/02_bam/8851598-Blood.s.md.ir.br.bam /home/users/team_projects/uveal_melanoma/02_bam/8834201-Blood.s.md.ir.br.bam /home/users/team_projects/uveal_melanoma/02_bam/8830513-Blood.s.md.ir.br.bam /home/users/team_projects/uveal_melanoma/02_bam/8139137-Blood.s.md.ir.br.bam /home/users/team_projects/uveal_melanoma/02_bam/8837700-Blood.s.md.ir.br.bam /home/users/team_projects/uveal_melanoma/02_bam/3638927-Blood.s.md.ir.br.bam /home/users/team_projects/uveal_melanoma/02_bam/8837803-Blood.s.md.ir.br.bam /home/users/team_projects/uveal_melanoma/02_bam/8859677-Blood.s.md.ir.br.bam /home/users/team_projects/uveal_melanoma/02_bam/8856619-Blood.s.md.ir.br.bam /home/users/team_projects/uveal_melanoma/02_bam/8852025-Blood.s.md.ir.br.bam /home/users/team_projects/uveal_melanoma/02_bam/10024819-Blood.s.md.ir.br.bam /home/users/team_projects/uveal_melanoma/02_bam/10007708-Blood.s.md.ir.br.bam /home/users/team_projects/uveal_melanoma/02_bam/7900560-Blood.s.md.ir.br.bam /home/users/team_projects/uveal_melanoma/02_bam/9151344-Blood.s.md.ir.br.bam /home/users/team_projects/uveal_melanoma/02_bam/10017731-Blood.s.md.ir.br.bam /home/users/team_projects/uveal_melanoma/02_bam/10038025-Blood.s.md.ir.br.bam /home/users/team_projects/uveal_melanoma/02_bam/10032815-Blood.s.md.ir.br.bam -o korean36.36s.q0Q0.chr20.mpileup) 2> korean36.mpu.chr20.out && mv korean36.mpu.chr20.out korean36.mpu.chr20.success