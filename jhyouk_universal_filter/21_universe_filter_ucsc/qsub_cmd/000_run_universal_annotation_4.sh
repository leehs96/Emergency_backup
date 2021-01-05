#!/bin/bash
#PBS -l nodes=1:ppn=1
#PBS -j oe
#PBS -o /home/users/jhyouk/81_filter_test_LADC/11_universe_filter/qsub_sdout/000_run_universal_annotation_4.sh.sdout
cd /home/users/jhyouk/81_filter_test_LADC/11_universe_filter
#./000_universal_annotation_filter.sh TCGA-05-4397.tumor indel hg19 /home/users/team_projects/LungAdeno_WGS_sypark/02_BAM/TCGA-05-4397.tumor.bam /home/users/team_projects/LungAdeno_WGS_sypark/02_BAM/TCGA-05-4397.normal.bam ../10_PanalOfNormal/pcawg3.3s.q0Q0.chr1.mpileup.indel.edit.gz