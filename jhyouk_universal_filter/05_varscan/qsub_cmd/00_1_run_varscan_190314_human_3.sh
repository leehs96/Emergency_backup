#!/bin/bash
#PBS -l nodes=1:ppn=1
#PBS -j oe
#PBS -o /home/users/jhyouk/81_filter_test_LADC/05_varscan/qsub_sdout/00_1_run_varscan_190314_human_3.sh.sdout
cd /home/users/jhyouk/81_filter_test_LADC/05_varscan
sh 00_human_b37_varscan.sh ../03_mpileup TCGA-55-6972.tumor TCGA-55-6972.normal