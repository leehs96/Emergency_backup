#!/bin/bash
#PBS -l nodes=1:ppn=1
#PBS -j oe
#PBS -o /home/users/jhyouk/81_filter_test_LADC/10_PanalOfNormal/qsub_sdout/02_1_run_edit_snp_36_190917_15.sh.sdout
cd /home/users/jhyouk/81_filter_test_LADC/10_PanalOfNormal
python 02_normal_panel_snv_edit.py korean36.36s.q0Q0.chr15.mpileup.snv.gz 36