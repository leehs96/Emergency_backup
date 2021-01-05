#!/bin/bash
#PBS -l nodes=1:ppn=1
#PBS -j oe
#PBS -o /home/users/jhyouk/81_filter_test_LADC/10_PanalOfNormal/qsub_sdout/temp_1.sh.sdout
cd /home/users/jhyouk/81_filter_test_LADC/10_PanalOfNormal
gzip -v -t korean19.19s.q0Q0.chr1.mpileup.snv.gz &> chr1.snv.gzip.check.out