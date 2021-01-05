#!/bin/bash
#PBS -l nodes=1:ppn=1
#PBS -j oe
#PBS -o /home/users/jhyouk/81_filter_test_LADC/10_PanalOfNormal/qsub_sdout/02_1_run_edit_indel_3.sh.sdout
cd /home/users/jhyouk/81_filter_test_LADC/10_PanalOfNormal
python 02_normal_panel_indel_edit_181202.py pcawg3.3s.q0Q0.chr3.mpileup.indel.gz 3