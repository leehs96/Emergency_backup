#!/bin/bash
#PBS -l nodes=1:ppn=1
#PBS -j oe
#PBS -o /home/users/jhyouk/81_filter_test_LADC/10_PanalOfNormal/qsub_sdout/run.pcawg3.mpu_14.sh.sdout
cd /home/users/jhyouk/81_filter_test_LADC/10_PanalOfNormal
(samtools mpileup -AB -d 3000 -q 0 -Q 0 -f /home/users/jhyouk/99_reference/human/GRCh37/human_g1k_v37.fasta -r 17 /home/users/team_projects/LungAdeno_WGS_sypark/02_BAM/TCGA-05-4397.normal.bam /home/users/team_projects/LungAdeno_WGS_sypark/02_BAM/TCGA-50-6591.normal.bam /home/users/team_projects/LungAdeno_WGS_sypark/02_BAM/TCGA-55-6972.normal.bam -o pcawg3.3s.q0Q0.chr17.mpileup) 2> pcawg3.mpu.chr17.out && mv pcawg3.mpu.chr17.out pcawg3.mpu.chr17.success