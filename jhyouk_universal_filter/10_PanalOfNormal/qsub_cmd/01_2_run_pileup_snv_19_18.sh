#!/bin/bash
#PBS -l nodes=1:ppn=1
#PBS -j oe
#PBS -o /home/users/jhyouk/81_filter_test_LADC/10_PanalOfNormal/qsub_sdout/01_2_run_pileup_snv_19_18.sh.sdout
cd /home/users/jhyouk/81_filter_test_LADC/10_PanalOfNormal
python 01_pileup_calling_snv.py korean19.19s.q0Q0.chr18.mpileup 19