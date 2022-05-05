################################
# PRS-CSx: population assignment
################################


# - Sex check (“Submitted.Gender” vs. “Inferred.Gender”)
# - Missingness/Heterozygosity (“het.missing.outliers”)
# - Sex chromosome aneuploidy (“putative.sex.chromosome.aneuploidy”)
# - Relatedness filter: see extr_unrel_inds.sh

# module load R/3.4.0
# source ~/R-3.4.0-ownlib.bash
# R

library(data.table)
library(dplyr)


# Read in data
eid <- fread("/data/tge/Tian/UKBB_full/genomics/cal_raw/ukb_cal_chr22_v2.fam", h=F)
eid$index <- 1:nrow(eid) # 488377
ukbb.sqc <- fread("/data/tge/Tian/UKBB_full/genomics/qc/ukb_sqc_v2.txt", h=F) # 488377
ukbb.sqc <- cbind(eid$V1, ukbb.sqc)
names(ukbb.sqc)[1] <- "eid"

ukbb.pred.pop <- fread("/data/js95/yfeng/projects/prscsx_pca/results/ukbb_ref.projected.PC.w.pop.pred.tsv", h=T)
ukbb.sqc.in.neale <- fread("/data/js95/yfeng/projects/ukbb_aoo/genetic_data/ukb32568.gwas_samples.txt", h=F)
ukbb.rel <- fread("/data/js95/yfeng/projects/prscsx_pca/results/ukbb_rel_sample.remove.list", h=F)

# "Submitted.Gender"
table(ukbb.sqc$V10)
#      F      M
# 264772 223605
# "Inferred.Gender"
table(ukbb.sqc$V11)
#      F      M
# 264864 223513
table(ukbb.sqc$V10==ukbb.sqc$V11)
 # FALSE   TRUE
 #   378 487999

# “het.missing.outliers”
table(ukbb.sqc$V19)
#      0      1
# 487409    968

# “putative.sex.chromosome.aneuploidy”
table(ukbb.sqc$V20)
#      0      1
# 487725    652



# Keep samples that pass each QC filter
ukbb.sqc.clean <- ukbb.sqc %>% filter(V10 == V11 & V19 == 0 & V20 == 0) #sex check; het.missing.outliers; putative.sex.chromosome.aneuploidy; 486565
ukbb.sqc.clean <- ukbb.sqc.clean %>% filter(!eid %in% ukbb.rel$V1) #411913
ukbb.sqc.clean <- merge(ukbb.sqc.clean, ukbb.pred.pop[,c("IID","predicted_pop")], by.x="eid", by.y="IID")

# Check any eids that have negative values (consent withdrawn)
# sum(ukbb.sqc.clean$eid <0) #14
ukbb.sqc.clean <- ukbb.sqc.clean %>% filter(eid > 0) # --> 411899

# Check N in each pop after sample QC 
ukbb.sqc.clean.eur <- ukbb.sqc.clean %>% filter(predicted_pop=="EUR") #379691
ukbb.sqc.clean.afr <- ukbb.sqc.clean %>% filter(predicted_pop=="AFR") #7627
ukbb.sqc.clean.amr <- ukbb.sqc.clean %>% filter(predicted_pop=="AMR") #696
ukbb.sqc.clean.sas <- ukbb.sqc.clean %>% filter(predicted_pop=="SAS") #8524
ukbb.sqc.clean.eas <- ukbb.sqc.clean %>% filter(predicted_pop=="EAS") #2222
ukbb.sqc.clean.other <- ukbb.sqc.clean %>% filter(predicted_pop=="Other") #13139


# Intersect with Neale lab GWAS samples
ukbb.sqc.in.neale <- fread("/data/js95/yfeng/projects/ukbb_aoo/genetic_data/ukb32568.gwas_samples.txt", h=F)

# Intersection between each unrel predicted pop and Neale lab GWAS samples
table(ukbb.sqc.clean.eur$eid %in% ukbb.sqc.in.neale$V1)
 # FALSE   TRUE
 # 45941 333750
table(ukbb.sqc.clean.afr$eid %in% ukbb.sqc.in.neale$V1) # all F
table(ukbb.sqc.clean.amr$eid %in% ukbb.sqc.in.neale$V1) # all F
table(ukbb.sqc.clean.sas$eid %in% ukbb.sqc.in.neale$V1) # all F
table(ukbb.sqc.clean.eas$eid %in% ukbb.sqc.in.neale$V1) # all F
table(ukbb.sqc.clean.other$eid %in% ukbb.sqc.in.neale$V1) # 12993 F, 146 T


# Write out final pop labels
fwrite(ukbb.sqc.clean.eur[ukbb.sqc.clean.eur$eid %in% ukbb.sqc.in.neale$V1, c("eid","eid")], "/data/js95/yfeng/projects/prscsx_pca/results/ukbb_sqc_eur.in.neale2_predProb0.9.tsv", quote=F, col.names=F, row.names=F, sep="\t")
fwrite(ukbb.sqc.clean.eur[!ukbb.sqc.clean.eur$eid %in% ukbb.sqc.in.neale$V1, c("eid","eid")], "/data/js95/yfeng/projects/prscsx_pca/results/ukbb_sqc_eur.non.neale2_predProb0.9.tsv", quote=F, col.names=F, row.names=F, sep="\t")

fwrite(ukbb.sqc.clean.eur[,c("eid","eid")], "/data/js95/yfeng/projects/prscsx_pca/results/ukbb_sqc_eur2_predProb0.9.tsv", quote=F, col.names=F, row.names=F, sep="\t")
fwrite(ukbb.sqc.clean.afr[,c("eid","eid")], "/data/js95/yfeng/projects/prscsx_pca/results/ukbb_sqc_afr2_predProb0.9.tsv", quote=F, col.names=F, row.names=F, sep="\t")
fwrite(ukbb.sqc.clean.amr[,c("eid","eid")], "/data/js95/yfeng/projects/prscsx_pca/results/ukbb_sqc_amr2_predProb0.9.tsv", quote=F, col.names=F, row.names=F, sep="\t")
fwrite(ukbb.sqc.clean.sas[,c("eid","eid")], "/data/js95/yfeng/projects/prscsx_pca/results/ukbb_sqc_sas2_predProb0.9.tsv", quote=F, col.names=F, row.names=F, sep="\t")
fwrite(ukbb.sqc.clean.eas[,c("eid","eid")], "/data/js95/yfeng/projects/prscsx_pca/results/ukbb_sqc_eas2_predProb0.9.tsv", quote=F, col.names=F, row.names=F, sep="\t")

fwrite(ukbb.sqc.clean.other[,c("eid","eid")], "/data/js95/yfeng/projects/prscsx_pca/results/ukbb_sqc_other2_predProb0.9.tsv", quote=F, col.names=F, row.names=F, sep="\t")






