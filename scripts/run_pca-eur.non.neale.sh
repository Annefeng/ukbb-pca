#!/bin/bash
#BSUB -J run_pop_pca
#BSUB -o /data/js95/yfeng/projects/prscsx_pca/run_eur.non.neale_pca_91k.out
#BSUB -q big
#BSUB -W 24:00
#BSUB -R rusage[mem=10000]

wdir=/data/js95/yfeng/projects/prscsx_pca/results
PLINK=/data/js95/yfeng/software/plink
PLINK2=/data/js95/yfeng/software/plink2
fhighld_region=/data/js95/yfeng/projects/pbk_genomics_qc/misc/long_range_LD_intervals.txt
scr_find_atgc=/data/js95/yfeng/projects/pbk_genomics_qc/scripts/find_atgc_snps.py



# # LD pruning
# $PLINK \
# --bfile $wdir/ukb_cal_v2 \
# --keep /data/js95/yfeng/projects/ukbb_aoo/misc_data/ukb32568_20190412_eur_genomics.fam \
# --extract $wdir/ukb_cal_v2.bim.non-atgc.snplist \
# --autosome \
# --geno 0.02 \
# --maf 0.05 \
# --snps-only just-acgt \
# --exclude range $fhighld_region \
# --indep-pairwise 100 50 0.2 \
# --out $wdir/ukbb_eur.non.neale_predProb0.9.autsnp.geno02.maf05-ldpr


# Run PCA
$PLINK2 \
--bfile $wdir/ukb_cal_v2 \
--keep /data/js95/yfeng/projects/ukbb_aoo/misc_data/ukb32568_20190412_eur_genomics_nonNeale.fam \
--extract $wdir/ukbb_eur.non.neale_predProb0.9.autsnp.geno02.maf05-ldpr.prune.in \
--memory 10000 \
--pca 20 approx \
--out $wdir/ukbb_eur.non.neale_predProb0.9.autsnp.geno02.maf05.pca


