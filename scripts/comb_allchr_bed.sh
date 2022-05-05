#!/bin/bash
#BSUB -J make_allchr_bed
#BSUB -o /data/js95/yfeng/projects/prscsx_pca
#BSUB -q big
#BSUB -W 15:00
#BSUB -R rusage[mem=10000]

wdir=/data/js95/yfeng/projects/prscsx_pca/results
PLINK=/data/js95/yfeng/software/plink

cd $wdir
$PLINK \
--merge-list ukb_cal_chr_merge.txt \
--make-bed \
--out ukb_cal_v2
