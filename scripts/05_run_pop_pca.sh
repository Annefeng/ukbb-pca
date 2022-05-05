#!/bin/bash
#BSUB -J run_pop_pca[1-5]
#BSUB -o /data/js95/yfeng/projects/prscsx_pca/run_pop_pca-%J-%I.out
#BSUB -q big
#BSUB -W 3:00
#BSUB -R rusage[mem=5000]


if [ ${LSB_JOBINDEX} == 1 ]; then 
    pop=eur.non.neale
elif [ ${LSB_JOBINDEX} == 2 ]; then 
    pop=afr
elif [ ${LSB_JOBINDEX} == 3 ]; then 
    pop=amr
elif [ ${LSB_JOBINDEX} == 4 ]; then 
    pop=eas
elif [ ${LSB_JOBINDEX} == 5 ]; then 
    pop=sas
elif [ ${LSB_JOBINDEX} == 6 ]; then 
    pop=eur.non.neale2
fi

wdir=/data/js95/yfeng/projects/prscsx_pca/results
PLINK=/data/js95/yfeng/software/plink
PLINK2=/data/js95/yfeng/software/plink2
fhighld_region=/data/js95/yfeng/projects/pbk_genomics_qc/misc/long_range_LD_intervals.txt
scr_find_atgc=/data/js95/yfeng/projects/pbk_genomics_qc/scripts/find_atgc_snps.py


# Run pop-specific PCA (based on genotype data)
# data is stored by chromosome
# Combine all chrs
# pop-speciifc LD pruning
# pop-specific pca


# # Find strand ambiguous SNPs
# python $scr_find_atgc $wdir/ukb_cal_v2.bim > $wdir/ukb_cal_v2.bim.atgc.snplist


# # Write a list of non-strand ambiguous SNPs to keep
# awk 'NR==FNR{a[$1];next} !($2 in a) {print $2}' $wdir/ukb_cal_v2.bim.atgc.snplist $wdir/ukb_cal_v2.bim > $wdir/ukb_cal_v2.bim.non-atgc.snplist


# # LD pruning
# $PLINK \
# --bfile $wdir/ukb_cal_v2 \
# --keep ukbb_sqc_eur_predProb0.9-in-neale-no.tsv \
# --extract $wdir/ukb_cal_v2.bim.non-atgc.snplist \
# --autosome \
# --geno 0.02 \
# --maf 0.05 \
# --snps-only just-acgt \
# --exclude range $fhighld_region \
# --indep-pairwise 100 50 0.2 \
# --out $wdir/ukbb_sqc_eur.non.neale_predProb0.9.autsnp.geno02.maf05-ldpr
# 122576 SNPs left


# pop=eur.non.neale2 #when used my own relatedness filter (extr_unrel_inds.sh)
 
# for pop in eur.non.neale afr amr eas sas; do

    echo $pop

    if [ -f $wdir/ukbb_sqc_${pop}_predProb0.9.tsv ] && [ ! -f $wdir/ukbb_sqc_${pop}_predProb0.9.autsnp.geno02.maf05-ldpr.prune.in ]; then
        # LD pruning
        $PLINK \
        --bfile $wdir/ukb_cal_v2 \
        --keep $wdir/ukbb_sqc_${pop}_predProb0.9.tsv \
        --extract $wdir/ukb_cal_v2.bim.non-atgc.snplist \
        --autosome \
        --geno 0.02 \
        --maf 0.05 \
        --snps-only just-acgt \
        --exclude range $fhighld_region \
        --indep-pairwise 100 50 0.2 \
        --out $wdir/ukbb_sqc_${pop}_predProb0.9.autsnp.geno02.maf05-ldpr
    fi

    if [ -f $wdir/ukbb_sqc_${pop}_predProb0.9.autsnp.geno02.maf05-ldpr.prune.in ]; then
        # Run PCA
        $PLINK2 \
        --bfile $wdir/ukb_cal_v2 \
        --keep $wdir/ukbb_sqc_${pop}_predProb0.9.tsv \
        --extract $wdir/ukbb_sqc_${pop}_predProb0.9.autsnp.geno02.maf05-ldpr.prune.in \
        --memory 6000 \
        --pca 20 approx \
        --out $wdir/ukbb_sqc_${pop}_predProb0.9.autsnp.geno02.maf05.pca
                # --pca var-wts 20 \
    fi

# done


# r2 < 0.2
 # 167508 ukbb_sqc_afr_predProb0.9.autsnp.geno02.maf05-ldpr.prune.in
 # 120312 ukbb_sqc_amr_predProb0.9.autsnp.geno02.maf05-ldpr.prune.in
 #  90258 ukbb_sqc_eas_predProb0.9.autsnp.geno02.maf05-ldpr.prune.in
 # 122576 ukbb_sqc_eur.non.neale_predProb0.9.autsnp.geno02.maf05-ldpr.prune.in
 # 124428 ukbb_sqc_sas_predProb0.9.autsnp.geno02.maf05-ldpr.prune.in


# r2 < 0.1
 # 111785 ukbb_sqc_afr_predProb0.9.autsnp.geno02.maf05-ldpr.prune.in
 #  66111 ukbb_sqc_amr_predProb0.9.autsnp.geno02.maf05-ldpr.prune.in
 #  54213 ukbb_sqc_eas_predProb0.9.autsnp.geno02.maf05-ldpr.prune.in
 #  69954 ukbb_sqc_eur.non.neale_predProb0.9.autsnp.geno02.maf05-ldpr.prune.in
 #  72608 ukbb_sqc_sas_predProb0.9.autsnp.geno02.maf05-ldpr.prune.in



###### Run pca with 1kg

onekg=/data/js95/yfeng/utility/1kG/ALL.1KG_phase3.20130502.genotypes.maf005
awk '$1=="CEU" || $1=="TSI" || $1=="FIN" || $1=="GBR" || $1=="IBS"' $onekg.fam | awk '{print $1"\t"$2}' > $wdir/onekg_eur.fam

sort $wdir/ukbb_sqc_eur.non.neale_predProb0.9.autsnp.geno02.maf05-ldpr.prune.in > $wdir/ukbb_eur.bim.tmp
awk '{print $2}' $onekg.bim | sort > $wdir/onekg.bim.tmp
join $wdir/ukbb_eur.bim.tmp $wdir/onekg.bim.tmp > $wdir/ukbb_eur.prune.in.1kg.snps #120986
rm *.bim.tmp

$PLINK \
--bfile $onekg \
--keep-allele-order \
--keep $wdir/onekg_eur.fam \
--extract $wdir/ukbb_eur.prune.in.1kg.snps \
--make-bed \
--out onekg_eur

$PLINK \
--bfile $wdir/ukb_cal_v2 \
--keep-allele-order \
--keep $wdir/ukbb_sqc_eur.non.neale_predProb0.9.tsv \
--extract $wdir/ukbb_eur.prune.in.1kg.snps \
--make-bed \
--out ukbb_eur

$PLINK --bfile $wdir/ukbb_eur \
--keep-allele-order \
--bmerge $wdir/onekg_eur \
--make-bed \
--out $wdir/ukbb.1kg_eur


$PLINK --bfile $wdir/ukbb.1kg_eur \
--pca 20 header tabs \
--out $wdir/ukbb.1kg_eur.pca


