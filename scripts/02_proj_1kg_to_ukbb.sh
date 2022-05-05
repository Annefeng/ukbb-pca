
################################
# PRS-CSx: population assignment
################################

wdir=/data/js95/yfeng/projects/prscsx_pca/results
PLINK=/data/js95/yfeng/software/plink
PLINK2=/data/js95/yfeng/software/plink2



# Project PCs onto study sampels: see https://www.cog-genomics.org/plink/2.0/score#pca_project
# ukbb_gen=/data/tge/Tian/UKBB_full/genomics/cal_raw

### [Use the chr-combined ukbb genotype data]
cd $wdir
ls /data/tge/Tian/UKBB_full/genomics/cal_raw/ukb_cal_chr*_v2.bed > ukb_cal_chr_merge-bed.tmp
ls /data/tge/Tian/UKBB_full/genomics/cal_raw/ukb_cal_chr*_v2.bim > ukb_cal_chr_merge-bim.tmp
ls /data/tge/Tian/UKBB_full/genomics/cal_raw/ukb_cal_chr*_v2.fam > ukb_cal_chr_merge-fam.tmp
paste ukb_cal_chr_merge-bed.tmp ukb_cal_chr_merge-bim.tmp ukb_cal_chr_merge-fam.tmp > ukb_cal_chr_merge.txt
rm *merge*tmp

$PLINK \
--merge-list ukb_cal_chr_merge.txt \
--make-bed \
--out ukb_cal_v2

$PLINK2 \
--bfile ukb_cal_v2 \
--read-freq $wdir/ALL.1KG_phase3.20130502.autosome.snp.geno02.maf05.pca.acount \
--score $wdir/ALL.1KG_phase3.20130502.autosome.snp.geno02.maf05.pca.eigenvec.var 2 3 header-read no-mean-imputation variance-standardize list-variants \
--score-col-nums 5-24 \
--out $wdir/ukb_cal_proj_pcs-all-chr

# after double checking with the pca results based on chr-specifiic files:
# rm ukb_cal_v2.*

### [Use the original by-chr genoypte data]
for chr in {22..1}; do
    echo chr$chr
    $PLINK2 \
    --bfile /data/tge/Tian/UKBB_full/genomics/cal_raw/ukb_cal_chr${chr}_v2 \
    --read-freq $wdir/ALL.1KG_phase3.20130502.autosome.snp.geno02.maf05.pca.acount \
    --score $wdir/ALL.1KG_phase3.20130502.autosome.snp.geno02.maf05.pca.eigenvec.var 2 3 header-read no-mean-imputation variance-standardize list-variants cols=+scoresums \
    --score-col-nums 5-24 \
    --out $wdir/ukb_cal_proj_pcs_chr${chr}
done

### test: project 1kg sampels into 1kg pc space
onekg_snp=/data/js95/yfeng/utility/1kG/ALL.1KG_phase3.20130502.genotypes.maf005 #with MAF >= 0.5%
$PLINK2 \
--bfile $onekg_snp \
--read-freq $wdir/ALL.1KG_phase3.20130502.autosome.snp.geno02.maf05.pca.acount \
--score $wdir/ALL.1KG_phase3.20130502.autosome.snp.geno02.maf05.pca.eigenvec.var 2 3 header-read no-mean-imputation variance-standardize list-variants \
--score-col-nums 5-24 \
--out $wdir/test_onekg_proj_pcs



# 'variance-standardize' linearly transforms each variant's dosage vector to have mean zero, variance 1.
# "Note that these PCs will be scaled a bit differently from ref_data.eigenvec; you need to multiply or divide the PCs by a multiple of sqrt(eigenvalue) to put them on the same scale."
# see: https://groups.google.com/forum/#!topic/plink2-users/W6DL5-hs_Q4
# https://groups.google.com/forum/#!searchin/plink2-users/pca$20projection%7Csort:date/plink2-users/ct8owo91Syc/pLGHtRwkl2wJ

# Sum over all PC score across 22 chrs for each individual (for each PC: col 6-25 in each sscore by chr file)

# wc -l ALL.1KG_phase3.20130502.autosome.snp.geno02.maf05.pca.eigenvec.var
# same as the total of
# wc -l ukb_cal_proj_pcs_*.sscore.vars

# nsnp=$(tail -n+2 ALL.1KG_phase3.20130502.autosome.snp.geno02.maf05.pca.eigenvec.var | wc -l)
# for pc in 1; do 
for pc in {2..20}; do 
    echo "PC"$pc
    # sum over all chr-specific prs into one
    # Once the pheno for each chr is generated, sum them up...
    pc_col=$(expr $pc + 25)

    for x in {1..22}; do 
        awk -v pc=$pc_col 'NR>1{print $pc}' ukb_cal_proj_pcs_chr${x}.sscore > ukb_cal_proj_pc${pc}_chr${x}.sscore.tmp
    done
    
    paste -d " " ukb_cal_proj_pc${pc}_chr{1..22}.sscore.tmp > ukb_cal_proj_pc${pc}.sscore.tmp
    awk '{for(i=1;i<=NF;i++) j+=$i; printf("%.10f\n",j); j=0}' ukb_cal_proj_pc${pc}.sscore.tmp > ukb_cal_proj_pc${pc}.sscore


    # awk -v x=$nsnp '{print $1/$x}' ukb_cal_proj_pc${pc}.sscore #using this in awk gives "fatal: division by zero attempted" error
    # awk -v x=$nsnp '{print $1/(149501*2)}' ukb_cal_proj_pc${pc}.sscore > ukb_cal_proj_pc${pc}.avg.sscore

    rm ukb_cal_proj_pc${pc}*.tmp

done
# note that this number isn't right: 149501*2; b/c there are some missing values; need to use NMISS_ALLELE_CT in each chr-specifis sscore-> see below (can compare with: combine all chr bed files into one and run PCA, to avoid erros)

for x in {1..22}; do
    awk 'NR>1{print $4}' ukb_cal_proj_pcs_chr${x}.sscore > ukb_cal_proj_pcs_chr${x}.sscore.nmiss.tmp
done
paste -d " " ukb_cal_proj_pcs_chr{1..22}.sscore.nmiss.tmp > ukb_cal_proj_pcs.sscore.nmiss.tmp
awk '{for(i=1;i<=NF;i++) j+=$i; printf("%.10f\n",j); j=0}' ukb_cal_proj_pcs.sscore.nmiss.tmp > ukb_cal_proj_pcs.sscore.nmiss
rm ukb_cal_proj_pcs*.sscore.nmiss.tmp


paste -d " " ukb_cal_proj_pc{1..20}.sscore > ukb_cal_proj_pcs.sscore
# paste -d " " ukb_cal_proj_pc*.avg.sscore > ukb_cal_proj_pcs.avg.sscore
# rm ukb_cal_proj_pc{1..20}.sscore ukb_cal_proj_pc{1..20}.avg.sscore


# Convert projected PCs to the original PC scale:
# 1. divdie by total allele count: Nmiss alleles
# 2. scaling:
# see https://groups.google.com/forum/#!topic/plink2-users/W6DL5-hs_Q4
cd $wdir
module load R/3.4.0
source ~/R-3.4.0-ownlib.bash
R
###
library(data.table)


setwd("/data/js95/yfeng/projects/prscsx_pca/results")

### test: 1kg data as the target sample
prj_pca_raw <- fread("test_onekg_proj_pcs.sscore", h=T)
eid <- prj_pca_raw$IID
prj_pca_raw <- prj_pca_raw[,-c(1:4)]

eigenval <- fread("ALL.1KG_phase3.20130502.autosome.snp.geno02.maf05.pca.eigenval", h=F)
eigenval <- t(as.data.frame(eigenval))

prj_pca_raw <- as.data.frame(prj_pca_raw)
prj_pca_adj <- data.frame(matrix(ncol=ncol(prj_pca_raw), nrow=nrow(prj_pca_raw), NA))
for(i in 1:20){
    prj_pca_adj[,i] <- prj_pca_raw[,i] / (-sqrt(eigenval[i])/2)
}
prj_pca_adj <- cbind(eid, prj_pca_adj)
names(prj_pca_adj)[-1] <- paste0("PC",1:20)
head(prj_pca_adj,3)
true_pca <- fread("ALL.1KG_phase3.20130502.autosome.snp.geno02.maf05.pca.eigenvec", h=T)
head(true_pca,3)
summary(prj_pca_adj$PC1 - true_pca$PC1)
#       Min.    1st Qu.     Median       Mean    3rd Qu.       Max.
# -1.398e-06 -8.781e-07  3.986e-07  1.328e-08  4.636e-07  6.685e-07
summary(prj_pca_adj$PC2 - true_pca$PC2)
#       Min.    1st Qu.     Median       Mean    3rd Qu.       Max.
# -4.770e-07 -2.085e-07 -9.562e-08 -4.516e-09  2.480e-07  3.778e-07



### [Use the chr-combined ukbb genotype data]
prj_pca_raw <- fread("ukb_cal_proj_pcs-all-chr.sscore", h=T)
eid <- prj_pca_raw$IID
prj_pca_raw <- prj_pca_raw[,-c(1:4)]
eigenval <- fread("ALL.1KG_phase3.20130502.autosome.snp.geno02.maf05.pca.eigenval", h=F)
eigenval <- t(as.data.frame(eigenval))

prj_pca_raw <- as.data.frame(prj_pca_raw)
prj_pca_adj <- data.frame(matrix(ncol=ncol(prj_pca_raw), nrow=nrow(prj_pca_raw), NA))
for(i in 1:20){
    prj_pca_adj[,i] <- prj_pca_raw[,i] / (-sqrt(eigenval[i])/2)
}
prj_pca_adj <- cbind(eid, prj_pca_adj)
names(prj_pca_adj)[-1] <- paste0("PC",1:20)
prj_pca_adj[prj_pca_adj$eid == "1757503", ]
# > prj_pca_adj[prj_pca_adj$eid == "1757503", ]
#     eid          PC1        PC2        PC3         PC4        PC5
# 1757503 -0.009925168 0.02716779 0.01012664 -0.01756293 0.00133453
#          PC6         PC7           PC8          PC9        PC10
# -0.004228557 0.001718196 -0.0003338972 -0.001892617 0.003108308
#        PC11        PC12        PC13         PC14        PC15         PC16
# -0.01475794 0.005023839 0.008419481 0.0004543103 0.004952547 -0.003349791
#        PC17         PC18         PC19         PC20
# -0.01407717 -0.004827319 0.0005201522 -0.001012727

# names(prj_pca_adj)[1] <- "IID"
fwrite(prj_pca_adj, "ukb_cal_proj_pcs-all-chr.avg.scaled.sscore", quote=F, col.names=T, row.names=F, sep="\t")


prj_pca_adj <- prj_pca_adj[order(prj_pca_adj$eid), ]

prj_pca_adj0 <- fread("ukb_cal_proj_pcs.avg.scaled.sscore", h=T)
prj_pca_adj0 <- prj_pca_adj0[order(prj_pca_adj0$eid), ]

prj_pca_adj <- as.data.frame(prj_pca_adj)
prj_pca_adj0 <- as.data.frame(prj_pca_adj0)
for( i in 2:21 ){
    print(paste0("PC", i-1))
    print(cor(prj_pca_adj0[,i], prj_pca_adj[,i]))
}
#--> all cor == 1!
# cor(prj_pca_adj0[,-1], prj_pca_adj[,-1]) #all cor == 1
#--> can use either file


### [Use the original by-chr genoypte data]
eid <- fread("ukb_cal_proj_pcs_chr22.sscore", h=T)
eid <- eid$IID
eigenval <- fread("ALL.1KG_phase3.20130502.autosome.snp.geno02.maf05.pca.eigenval", h=F)
prj_pca_raw <- fread("ukb_cal_proj_pcs.sscore", h=F)
prj_pca_nmiss <- fread("ukb_cal_proj_pcs.sscore.nmiss", h=F)

eigenval <- t(as.data.frame(eigenval))
prj_pca_raw <- as.data.frame(prj_pca_raw)
prj_pca_adj <- data.frame(matrix(ncol=ncol(prj_pca_raw), nrow=nrow(prj_pca_raw), NA))
for(i in 1:20){
    prj_pca_adj[,i] <- prj_pca_raw[,i] / (prj_pca_nmiss * (-sqrt(eigenval[i])/2))
    # prj_pca_adj[,i] <- prj_pca_raw[,i] / (-sqrt(eigenval[i])/2)    
}
prj_pca_adj <- cbind(eid, prj_pca_adj)
names(prj_pca_adj)[-1] <- paste0("PC",1:20)
fwrite(prj_pca_adj, "ukb_cal_proj_pcs.avg.scaled.sscore", quote=F, col.names=T, row.names=F, sep="\t")

#       IID          PC1        PC2         PC3         PC4           PC5
# 1 1757503 -0.009925163 0.02716771 0.010126632 -0.01756292  1.334525e-03
#            PC6           PC7           PC8           PC9        PC10
# 1 -0.004228558  1.718196e-03 -0.0003338993 -0.0018926126 0.003108315
###
