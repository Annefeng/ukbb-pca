
################################
# PRS-CSx: population assignment
################################

wdir=/data/js95/yfeng/projects/prscsx_pca/results
PLINK=/data/js95/yfeng/software/plink
PLINK2=/data/js95/yfeng/software/plink2
fhighld_region=/data/js95/yfeng/projects/pbk_genomics_qc/misc/long_range_LD_intervals.txt
scr_find_atgc=/data/js95/yfeng/projects/pbk_genomics_qc/scripts/find_atgc_snps.py

# ukbb_snp=/data/js95/yfeng/projects/ukbb_aoo/genetic_data/ukb31063.gwas_variants.autosomes.tsv
onekg_snp=/data/js95/yfeng/utility/1kG/ALL.1KG_phase3.20130502.genotypes.maf005 #with MAF >= 0.5%
# $PLINK \
# --bfile $onekg_snp \
# --freq \
# --out $onekg_snp-frq


# Intersect 1KG with UKBB genotyped SNPs
cat /data/tge/Tian/UKBB_full/genomics/cal_raw/ukb_cal_chr*_v2.bim > $wdir/ukb_cal_v2.bim
awk 'NR==FNR{a[$2];next} ($2 in a) {print $2}' $wdir/ukb_cal_v2.bim $onekg_snp.bim > $wdir/ALL.1KG_phase3.20130502.genotypes.maf005.in.ukbb.snplist

# awk 'NR==FNR{a[$2];next} ($2 in a) {print $2}' $ukbb_snp $onekg_snp.bim > $wdir/ALL.1KG_phase3.20130502.genotypes.maf005.in.neale.ukb.snplist
# both files: column2: rsID


# Find strand ambiguous SNPs
python $scr_find_atgc $onekg_snp.bim > $wdir/ALL.1KG_phase3.20130502.genotypes.maf005.atgc.snplist


# Write a list of non-strand ambiguous SNPs to keep
awk 'NR==FNR{a[$1];next} !($1 in a) {print $1}' $wdir/ALL.1KG_phase3.20130502.genotypes.maf005.atgc.snplist $wdir/ALL.1KG_phase3.20130502.genotypes.maf005.in.ukbb.snplist > $wdir/ALL.1KG_phase3.20130502.genotypes.maf005.in.ukbb.non-atgc.snplist


# LD pruning
$PLINK \
--bfile $onekg_snp \
--extract $wdir/ALL.1KG_phase3.20130502.genotypes.maf005.in.ukbb.non-atgc.snplist \
--autosome \
--geno 0.02 \
--maf 0.05 \
--snps-only just-acgt \
--exclude range $fhighld_region \
--indep-pairwise 100 50 0.2 \
--out $wdir/ALL.1KG_phase3.20130502.autosome.snp.geno02.maf05-ldpr
# 149501 variants left


# $PLINK \
# --bfile $onekg_snp \
# --extract $wdir/ALL.1KG_phase3.20130502.genotypes.maf005.in.ukbb.non-atgc.snplist \
# --autosome \
# --geno 0.01 \
# --maf 0.05 \
# --snps-only just-acgt \
# --exclude range $fhighld_region \
# --indep-pairwise 100 50 0.2 \
# --out $wdir/ALL.1KG_phase3.20130502.autosome.snp.geno01.maf05-ldpr


# PCA of 1kg samples based on pruned/LD-indep SNPs: save all loadings
# $PLINK \
# --bfile $onekg_snp \
# --extract $wdir/ALL.1KG_phase3.20130502.autosome.snp.geno02.maf05-ldpr.prune.in \
# --pca 20 header tabs var-wts \
# --out $wdir/ALL.1KG_phase3.20130502.autosome.snp.geno02.maf05.pca


# $PLINK2 --help pca
$PLINK2 \
--bfile $onekg_snp \
--extract $wdir/ALL.1KG_phase3.20130502.autosome.snp.geno02.maf05-ldpr.prune.in \
--freq counts \
--pca var-wts 20 \
--out $wdir/ALL.1KG_phase3.20130502.autosome.snp.geno02.maf05.pca


# $PLINK2 \
# --bfile $onekg_snp \
# --extract $wdir/ALL.1KG_phase3.20130502.autosome.snp.geno01.maf05.rsq0.1-ldpr.prune.in \
# --freq counts \
# --pca var-wts 20 \
# --out $wdir/ALL.1KG_phase3.20130502.autosome.snp.geno01.maf05.rsq0.1.pca


