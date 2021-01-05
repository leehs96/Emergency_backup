#!/bin/bash
#PBS -l nodes=1:ppn=1
#PBS -j oe
#PBS -o /home/users/jhyouk/81_filter_test_LADC/13_SV/qsub_sdout/00_2_run_PON_1.sh.sdout
cd /home/users/jhyouk/81_filter_test_LADC/13_SV
python 00_Making_PON_Delly.py 00_1_list_190907.txt &> PON_35.out