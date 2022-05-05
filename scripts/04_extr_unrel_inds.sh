
################################
# PRS-CSx: population assignment
################################


module load R/3.4.0
source ~/R-3.4.0-ownlib.bash
R

library(data.table)
library(dplyr)
# For EUR individuals: separate those in Neale GWAS vs. not

neale.ukbb.sample <- fread("/data/js95/yfeng/projects/prscsx_pca/misc/samples.both_sexes.tsv", h=T)
# 361194

eid <- fread("/data/tge/Tian/UKBB_full/genomics/cal_raw/ukb_cal_chr22_v2.fam", h=F)
eid$index <- 1:nrow(eid)
# 488377
tian.ukbb.sqc <- fread("/data/tge/Tian/UKBB_full/genomics/qc/ukb_sqc_v2.txt", h=F)
# 488377
tian.ukbb.sqc <- tian.ukbb.sqc[,c(5:6)]
names(tian.ukbb.sqc) <- names(neale.ukbb.sample)
tian.ukbb.sqc$index <- 1:nrow(tian.ukbb.sqc)

tian.ukbb.sqc.in.neale <- merge(tian.ukbb.sqc, neale.ukbb.sample, by=c("plate_name","well"))
nrow(tian.ukbb.sqc.in.neale) #361206 -> some duplicates

# > table(duplicated(neale.ukbb.sample))
#  FALSE   TRUE
# 361192      2

# > neale.ukbb.sample[duplicated(neale.ukbb.sample),]
#       plate_name well
# 1: SMP4_0014502A  E05
# 2: SMP4_0014640A  H04
# neale.ukbb.sample %>% filter(plate_name=="SMP4_0014502A" & well=="E05")

# tian.ukbb.sqc.in.neale[duplicated(tian.ukbb.sqc.in.neale$index),]
#       plate_name well  index
# 1: SMP4_0014502A  E05 221265
# 2: SMP4_0014502A  E05 240651
# 3: SMP4_0014640A  H04 228769
# 4: SMP4_0014640A  H04 435875

# > tian.ukbb.sqc.in.neale %>% filter(plate_name=="SMP4_0014502A" & well=="E05")
#      plate_name well  index
# 1 SMP4_0014502A  E05 221265
# 2 SMP4_0014502A  E05 221265
# 3 SMP4_0014502A  E05 240651
# 4 SMP4_0014502A  E05 240651
# > tian.ukbb.sqc.in.neale %>% filter(plate_name=="SMP4_0014640A" & well=="H04")
#      plate_name well  index
# 1 SMP4_0014640A  H04 228769
# 2 SMP4_0014640A  H04 228769
# 3 SMP4_0014640A  H04 435875
# 4 SMP4_0014640A  H04 435875

# tian.ukbb.sqc.in.neale %>% filter(index=="221265")
#      plate_name well  index
# 1 SMP4_0014502A  E05 221265
# 2 SMP4_0014502A  E05 221265


# tian.ukbb.sqc.in.neale %>% distinct() %>% dim() #361202


# > tian.ukbb.sqc.in.neale[duplicated(tian.ukbb.sqc.in.neale[,-3])]
#        plate_name well  index
#  1: SMP4_0012383A  C09 245682
#  2: SMP4_0013746A  H09 240254
#  3: SMP4_0014502A  A08 240650
#  4: SMP4_0014502A  E05 221265
#  5: SMP4_0014502A  E05 240651
#  6: SMP4_0014502A  E05 240651
#  7: SMP4_0014503A  F01 240652
#  8: SMP4_0014640A  H04 228769
#  9: SMP4_0014640A  H04 435875
# 10: SMP4_0014640A  H04 435875
# 11: SMP4_0014641A  B04 240775
# 12: SMP4_0014641A  C05 240776
# 13: SMP4_0016202A  B01 241453
# 14: SMP4_0016202A  C01 241454

# tian.ukbb.sqc.in.neale %>% filter(plate_name=="SMP4_0012383A" & well=="C09")
#      plate_name well  index
# 1 SMP4_0012383A  C09 195880
# 2 SMP4_0012383A  C09 245682

dup_index <- tian.ukbb.sqc.in.neale[,count := .N,by = list(plate_name,well)][count>1, which=T]
tian.ukbb.sqc.in.neale[dup_index, ]
#        plate_name well  index count
#  1: SMP4_0012383A  C09 195880     2
#  2: SMP4_0012383A  C09 245682     2
#  3: SMP4_0013746A  H09 198626     2
#  4: SMP4_0013746A  H09 240254     2
#  5: SMP4_0014502A  A08 221222     2
#  6: SMP4_0014502A  A08 240650     2
#  7: SMP4_0014502A  E05 221265     4
#  8: SMP4_0014502A  E05 221265     4
#  9: SMP4_0014502A  E05 240651     4
# 10: SMP4_0014502A  E05 240651     4
# 11: SMP4_0014503A  F01 221356     2
# 12: SMP4_0014503A  F01 240652     2
# 13: SMP4_0014640A  H04 228769     4
# 14: SMP4_0014640A  H04 228769     4
# 15: SMP4_0014640A  H04 435875     4
# 16: SMP4_0014640A  H04 435875     4
# 17: SMP4_0014641A  B04 228792     2
# 18: SMP4_0014641A  B04 240775     2
# 19: SMP4_0014641A  C05 228805     2
# 20: SMP4_0014641A  C05 240776     2
# 21: SMP4_0016202A  B01 236985     2
# 22: SMP4_0016202A  B01 241453     2
# 23: SMP4_0016202A  C01 236997     2
# 24: SMP4_0016202A  C01 241454     2

tian.ukbb.sqc.in.neale <- tian.ukbb.sqc.in.neale %>% distinct() #361202
dup_index <- tian.ukbb.sqc.in.neale[,count := .N,by = list(plate_name,well)][count>1, which=T]
tian.ukbb.sqc.in.neale[dup_index, ]
#        plate_name well  index count
#  1: SMP4_0012383A  C09 195880     2
#  2: SMP4_0012383A  C09 245682     2
#  3: SMP4_0013746A  H09 198626     2
#  4: SMP4_0013746A  H09 240254     2
#  5: SMP4_0014502A  A08 221222     2
#  6: SMP4_0014502A  A08 240650     2
#  7: SMP4_0014502A  E05 221265     2
#  8: SMP4_0014502A  E05 240651     2
#  9: SMP4_0014503A  F01 221356     2
# 10: SMP4_0014503A  F01 240652     2
# 11: SMP4_0014640A  H04 228769     2
# 12: SMP4_0014640A  H04 435875     2
# 13: SMP4_0014641A  B04 228792     2
# 14: SMP4_0014641A  B04 240775     2
# 15: SMP4_0014641A  C05 228805     2
# 16: SMP4_0014641A  C05 240776     2
# 17: SMP4_0016202A  B01 236985     2
# 18: SMP4_0016202A  B01 241453     2
# 19: SMP4_0016202A  C01 236997     2
# 20: SMP4_0016202A  C01 241454     2

tian.ukbb.sqc.in.neale <- merge(tian.ukbb.sqc.in.neale, eid[,c(1,7)], by="index")
names(tian.ukbb.sqc.in.neale)[5] <- "eid"


# EUR not in Neale GWAS
pop_label <- fread("/data/js95/yfeng/projects/prscsx_pca/results/ukbb_ref.projected.PC.w.pop.pred.tsv", h=T)
pop_label.eur <- pop_label[pop_label$predicted_pop=="EUR", ]
pop_label.eur$in.neale <- ifelse(pop_label.eur$IID %in% tian.ukbb.sqc.in.neale$eid, "YES", "NO")
# > table(pop_label.eur$in.neale)
#     NO    YES
#  92302 361042


tian.ukbb.sqc.in.neale.2 <- fread("/data/js95/yfeng/projects/ukbb_aoo/genetic_data/ukb32568.gwas_samples.txt", h=F)
pop_label.eur$in.neale.2 <- ifelse(pop_label.eur$IID %in% tian.ukbb.sqc.in.neale.2$V1, "YES", "NO")
# > table(pop_label.eur$in.neale.2)

#     NO    YES
#  92309 361035

# pop_label$in.neale <- ifelse(pop_label$IID %in% tian.ukbb.sqc.in.neale.2$V1, "YES", "NO")
# > table(pop_label$in.neale)
#     NO    YES
# 127183 361194
# > table(pop_label$predicted_pop[pop_label$in.neale == "YES"])
#    EUR  Other
# 361035    159


##############
# Extract unrelated inds for other pop groups
tian.ukbb.sqc <- fread("/data/tge/Tian/UKBB_full/genomics/qc/ukb_sqc_v2.txt", h=F)
# N related (used in pca calc- see Data-Field 22020 and supp info in Bycroft et al Nature)
table(tian.ukbb.sqc$V25)
 #     0      1
 # 81158 407219


tian.ukbb.sqc <- cbind(eid$V1, tian.ukbb.sqc)
names(tian.ukbb.sqc)[1] <- "eid"
tian.ukbb.sqc <- merge(tian.ukbb.sqc, pop_label[,c("IID","predicted_pop")], by.x="eid", by.y="IID")

tian.ukbb.sqc.unrel <- tian.ukbb.sqc %>% filter(V25 == 1) #488377 --> 407219
head(tian.ukbb.sqc.unrel)

# Check any eids that have negative values (consent withdrawn)
sum(tian.ukbb.sqc.unrel$eid <0) #13
# tian.ukbb.sqc.unrel <- tian.ukbb.sqc.unrel %>% filter(eid > 0)

tian.ukbb.sqc.unrel.eur <- tian.ukbb.sqc.unrel %>% filter(predicted_pop=="EUR") #375498
tian.ukbb.sqc.unrel.afr <- tian.ukbb.sqc.unrel %>% filter(predicted_pop=="AFR") #7513
tian.ukbb.sqc.unrel.amr <- tian.ukbb.sqc.unrel %>% filter(predicted_pop=="AMR") #689
tian.ukbb.sqc.unrel.sas <- tian.ukbb.sqc.unrel %>% filter(predicted_pop=="SAS") #8420
tian.ukbb.sqc.unrel.eas <- tian.ukbb.sqc.unrel %>% filter(predicted_pop=="EAS") #2181
tian.ukbb.sqc.unrel.other <- tian.ukbb.sqc.unrel %>% filter(predicted_pop=="Other") #12918



# Find unrelated EUR not in Neale
tian.ukbb.sqc.in.neale <- fread("/data/js95/yfeng/projects/ukbb_aoo/genetic_data/ukb32568.gwas_samples.txt", h=F)

# Intersection between each unrel predicted pop and Neale lab GWAS samples
table(tian.ukbb.sqc.unrel.eur$eid %in% tian.ukbb.sqc.in.neale$V1)
 # FALSE   TRUE
 # 14463 361035
table(tian.ukbb.sqc.unrel.afr$eid %in% tian.ukbb.sqc.in.neale$V1) # all F
table(tian.ukbb.sqc.unrel.amr$eid %in% tian.ukbb.sqc.in.neale$V1) # all F
table(tian.ukbb.sqc.unrel.sas$eid %in% tian.ukbb.sqc.in.neale$V1) # all F
table(tian.ukbb.sqc.unrel.eas$eid %in% tian.ukbb.sqc.in.neale$V1) # all F
table(tian.ukbb.sqc.unrel.other$eid %in% tian.ukbb.sqc.in.neale$V1) # 12759 F, 159 T



#####
# Do this staring from the relatedness file provided
# Remove one from each pair of related individuals (w/o removing duplicated ones)
ukbb.rel <- fread("/data/tge/Tian/UKBB_full/genomics/qc/ukb32568_rel_s488363.dat", h=T)
# UKBB Resource 531: This file lists the pairs of individuals related up to the third degree in the data set. It is a plaintext file with space separated columns.
samples <- names(sort(table(unlist(ukbb.rel[,c('ID1', 'ID2')])), decreasing=TRUE)) #len: 147724
removed <- character()

for (row in 1:dim(ukbb.rel)[1]) {
  
  print(row)
  i_sample <- as.character(ukbb.rel[row, 'ID1'])
  j_sample <- as.character(ukbb.rel[row, 'ID2'])
  
  i_index <- match(i_sample, samples)
  j_index <- match(j_sample, samples)
  
  if (is.na(i_index) | is.na(j_index)) { next }
  
  if (i_index <= j_index) {
    removed <- c(removed, i_sample)
    samples <- samples[!samples == i_sample]
  } else {
    removed <- c(removed, j_sample)
    samples <- samples[!samples == j_sample]
  }
  
}
length(removed) #74788
length(samples) #72936


write.table(removed, "/data/js95/yfeng/projects/prscsx_pca/results/ukbb_rel_sample.remove.list", quote=F, col.names=F, row.names=F)

# removed <- fread("/data/js95/yfeng/projects/prscsx_pca/results/ukbb_rel_sample.remove.list", h=F)
# removed <- removed$V1

tian.ukbb.sqc <- cbind(eid$V1, tian.ukbb.sqc)
names(tian.ukbb.sqc)[1] <- "eid"
tian.ukbb.sqc <- merge(tian.ukbb.sqc, pop_label[,c("IID","predicted_pop")], by.x="eid", by.y="IID")

tian.ukbb.sqc.unrel <- tian.ukbb.sqc %>% filter(!eid %in% removed) #488377 --> 413589
head(tian.ukbb.sqc.unrel)

# Remove any eids that have negative values (consent withdrawn)
tian.ukbb.sqc.unrel <- tian.ukbb.sqc.unrel %>% filter(eid > 0) #413575

tian.ukbb.sqc.unrel.eur <- tian.ukbb.sqc.unrel %>% filter(predicted_pop=="EUR") #381198
tian.ukbb.sqc.unrel.afr <- tian.ukbb.sqc.unrel %>% filter(predicted_pop=="AFR") #7635
tian.ukbb.sqc.unrel.amr <- tian.ukbb.sqc.unrel %>% filter(predicted_pop=="AMR") #699
tian.ukbb.sqc.unrel.sas <- tian.ukbb.sqc.unrel %>% filter(predicted_pop=="SAS") #8588
tian.ukbb.sqc.unrel.eas <- tian.ukbb.sqc.unrel %>% filter(predicted_pop=="EAS") #2225
tian.ukbb.sqc.unrel.other <- tian.ukbb.sqc.unrel %>% filter(predicted_pop=="Other") #13230



# Find unrelated EUR not in Neale
tian.ukbb.sqc.in.neale <- fread("/data/js95/yfeng/projects/ukbb_aoo/genetic_data/ukb32568.gwas_samples.txt", h=F)
tian.ukbb.sqc.unrel.eur$in.neale <- ifelse(tian.ukbb.sqc.unrel.eur$eid %in% tian.ukbb.sqc.in.neale$V1, "YES", "NO")

# Intersection between each unrel predicted pop and Neale lab GWAS samples
# > table(tian.ukbb.sqc.unrel.eur$in.neale)
#     NO    YES
#  47448 333750
table(tian.ukbb.sqc.unrel.afr$eid %in% tian.ukbb.sqc.in.neale$V1) # all F
table(tian.ukbb.sqc.unrel.amr$eid %in% tian.ukbb.sqc.in.neale$V1) # all F
table(tian.ukbb.sqc.unrel.sas$eid %in% tian.ukbb.sqc.in.neale$V1) # all F
table(tian.ukbb.sqc.unrel.eas$eid %in% tian.ukbb.sqc.in.neale$V1) # all F
table(tian.ukbb.sqc.unrel.other$eid %in% tian.ukbb.sqc.in.neale$V1) # 13084 F, 146 T


# Intersection between each predicted pop (incl. rel and unrel) and Neale lab GWAS samples
table(pop_label$predicted_pop[pop_label$IID %in% tian.ukbb.sqc.in.neale$V1])
#    EUR  Other
# 361035    159


fwrite(tian.ukbb.sqc.unrel.eur[tian.ukbb.sqc.unrel.eur$in.neale=="YES", c(1,5)], "/data/js95/yfeng/projects/prscsx_pca/results/ukbb_unrel_eur.in.neale_predProb0.9.tsv", quote=F, col.names=T, row.names=F, sep="\t")

fwrite(tian.ukbb.sqc.unrel.eur[tian.ukbb.sqc.unrel.eur$in.neale=="NO", c(1,5)], "/data/js95/yfeng/projects/prscsx_pca/results/ukbb_unrel_eur.non.neale_predProb0.9.tsv", quote=F, col.names=T, row.names=F, sep="\t")


fwrite(tian.ukbb.sqc.unrel.eur[,c(1,5)], "/data/js95/yfeng/projects/prscsx_pca/results/ukbb_unrel_eur_predProb0.9.tsv", quote=F, col.names=T, row.names=F, sep="\t")
fwrite(tian.ukbb.sqc.unrel.afr[,c(1,5)], "/data/js95/yfeng/projects/prscsx_pca/results/ukbb_unrel_afr_predProb0.9.tsv", quote=F, col.names=T, row.names=F, sep="\t")
fwrite(tian.ukbb.sqc.unrel.amr[,c(1,5)], "/data/js95/yfeng/projects/prscsx_pca/results/ukbb_unrel_amr_predProb0.9.tsv", quote=F, col.names=T, row.names=F, sep="\t")
fwrite(tian.ukbb.sqc.unrel.sas[,c(1,5)], "/data/js95/yfeng/projects/prscsx_pca/results/ukbb_unrel_sas_predProb0.9.tsv", quote=F, col.names=T, row.names=F, sep="\t")
fwrite(tian.ukbb.sqc.unrel.eas[,c(1,5)], "/data/js95/yfeng/projects/prscsx_pca/results/ukbb_unrel_eas_predProb0.9.tsv", quote=F, col.names=T, row.names=F, sep="\t")

fwrite(tian.ukbb.sqc.unrel.other[,c(1,5)], "/data/js95/yfeng/projects/prscsx_pca/results/ukbb_unrel_other_predProb0.9.tsv", quote=F, col.names=T, row.names=F, sep="\t")


