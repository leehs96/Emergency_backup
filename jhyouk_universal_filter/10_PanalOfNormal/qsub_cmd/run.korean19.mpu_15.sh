#!/bin/bash
#PBS -l nodes=1:ppn=1
#PBS -j oe
#PBS -o /home/users/jhyouk/81_filter_test_LADC/10_PanalOfNormal/qsub_sdout/run.korean19.mpu_15.sh.sdout
cd /home/users/jhyouk/81_filter_test_LADC/10_PanalOfNormal
(samtools mpileup -AB -d 3000 -q 0 -Q 0 -f /home/users/jhyouk/99_reference/human/GRCh37/human_g1k_v37.fasta -r 3 /home/users/team_projects/BRCA_yonsei/02_bam/4B2.s.md.ir.br.bam /home/users/team_projects/BRCA_yonsei/02_bam/ON015.s.md.ir.br.bam /home/users/team_projects/BRCA_yonsei/02_bam/2B1.s.md.ir.br.bam /home/users/team_projects/uveal_melanoma/02_bam/1-17-Blood.s.md.ir.br.bam /home/users/team_projects/uveal_melanoma/02_bam/2-24-Blood.s.md.ir.br.bam /home/users/team_projects/uveal_melanoma/02_bam/3-29-Blood.s.md.ir.br.bam /home/users/team_projects/uveal_melanoma/02_bam/4-42-Blood.s.md.ir.br.bam /home/users/team_projects/uveal_melanoma/02_bam/5-68-Normal.s.md.ir.br.bam /home/users/team_projects/uveal_melanoma/02_bam/6-88-Blood.s.md.ir.br.bam /home/users/team_projects/uveal_melanoma/02_bam/8802875-Blood.s.md.ir.br.bam /home/users/team_projects/uveal_melanoma/02_bam/8769056-Blood.s.md.ir.br.bam /home/users/team_projects/uveal_melanoma/02_bam/8804915-Blood.s.md.ir.br.bam /home/users/team_projects/uveal_melanoma/02_bam/8808738-Blood.s.md.ir.br.bam /home/users/team_projects/uveal_melanoma/02_bam/8812875-Blood.s.md.ir.br.bam /home/users/team_projects/Radiation_signature/05_bam_human_sample_GRCh37/HBIR1_germline.s.md.ir.br.bam /home/users/team_projects/Radiation_signature/05_bam_human_sample_GRCh37/HBIR2_germ.s.md.ir.br.bam /home/users/team_projects/Radiation_signature/05_bam_human_sample_GRCh37/HCIR1.s.md.ir.br.bam /home/users/team_projects/Radiation_signature/05_bam_human_sample_GRCh37/hs_study5_germline_HCIR2.s.md.ir.br.bam /home/users/team_projects/Radiation_signature/05_bam_human_sample_GRCh37/hs_study5_germline_HCNor1.s.md.ir.br.bam -o korean19.19s.q0Q0.chr3.mpileup) 2> korean19.mpu.chr3.out && mv korean19.mpu.chr3.out korean19.mpu.chr3.success