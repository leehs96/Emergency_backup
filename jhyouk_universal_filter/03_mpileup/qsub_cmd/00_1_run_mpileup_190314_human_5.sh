#!/bin/bash
#PBS -l nodes=1:ppn=1
#PBS -j oe
#PBS -o /home/users/jhyouk/81_filter_test_LADC/03_mpileup/qsub_sdout/00_1_run_mpileup_190314_human_5.sh.sdout
cd /home/users/jhyouk/81_filter_test_LADC/03_mpileup
sh 00_human_b37_script_mpileup.sh /home/users/team_projects/LungAdeno_WGS_sypark/02_BAM/TCGA-55-6972.tumor.bam TCGA-55-6972.tumor