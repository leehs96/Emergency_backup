#!/bin/bash
#PBS -l nodes=1:ppn=5
#PBS -j oe
#PBS -o /home/users/jhyouk/81_filter_test_LADC/01_delly/qsub_sdout/01_run_Delly_190314_human_3.sh.sdout
cd /home/users/jhyouk/81_filter_test_LADC/01_delly
sh 00_human_b37_Delly.sh /home/users/team_projects/LungAdeno_WGS_sypark/02_BAM/TCGA-55-6972.tumor.bam /home/users/team_projects/LungAdeno_WGS_sypark/02_BAM/TCGA-55-6972.normal.bam TCGA-55-6972.tumor