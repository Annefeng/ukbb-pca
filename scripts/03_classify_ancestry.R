
################################
# PRS-CSx: population assignment
################################

# module load R/3.4.0
# source ~/R-3.4.0-ownlib.bash
# R

setwd("/data/js95/yfeng/projects/prscsx_pca/results")

library(data.table)
library(randomForest)
library(tidyverse)
library(ggsci)
library(RColorBrewer)

pca_ukb <- fread("ukb_cal_proj_pcs.avg.scaled.sscore", h=T)
pca_ukb$FID <- "UKBB"
pca_ukb <- pca_ukb[, c(22,1:21)]
names(pca_ukb)[2] <- "IID"
pca_1kg <- fread("ALL.1KG_phase3.20130502.autosome.snp.geno02.maf05.pca.eigenvec", h=T)
names(pca_1kg)[1] <- "FID"

pca_comb <- rbind(pca_ukb, pca_1kg)


pca_comb$superpop <- "UKBB"
pca_comb$superpop <- ifelse(pca_comb$FID=="CHB" | pca_comb$FID=="JPT" | pca_comb$FID=="CHS" | pca_comb$FID=="CDX" | pca_comb$FID=="KHV", "EAS", pca_comb$superpop)
pca_comb$superpop <- ifelse(pca_comb$FID=="CEU" | pca_comb$FID=="TSI" | pca_comb$FID=="FIN" | pca_comb$FID=="GBR" | pca_comb$FID=="IBS", "EUR", pca_comb$superpop)
pca_comb$superpop <- ifelse(pca_comb$FID=="YRI" | pca_comb$FID=="LWK" | pca_comb$FID=="GWD" | pca_comb$FID=="MSL" | pca_comb$FID=="ESN" | pca_comb$FID=="ASW" | pca_comb$FID=="ACB", "AFR", pca_comb$superpop)
pca_comb$superpop <- ifelse(pca_comb$FID=="MXL" | pca_comb$FID=="PUR" | pca_comb$FID=="CLM" | pca_comb$FID=="PEL", "AMR", pca_comb$superpop)
pca_comb$superpop <- ifelse(pca_comb$FID=="GIH" | pca_comb$FID=="PJL" | pca_comb$FID=="BEB" | pca_comb$FID=="STU" | pca_comb$FID=="ITU", "SAS", pca_comb$superpop)

pca_comb$superpop <- factor(pca_comb$superpop,
                              levels=c("UKBB","AFR","AMR","EAS","EUR","SAS"))


# PC plots: projected ukbb pcs + 1kg pcs, colored by 1kg super population
for(i in 1:10){
    p = ggplot(pca_comb, aes_string(x=paste0('PC',i), y=paste0('PC',i+1))) +
    geom_point(aes(color=superpop), alpha=0.57, size=0.97) +
    scale_color_manual(values = c("grey25",pal_d3("category20")(20)[c(1:5)])) + 
    theme_bw() +
    guides(color = guide_legend(override.aes = list(size=2))) +
    labs(x=paste0('PC',i), y=paste0('PC',i+1), title="UKBB projected PCs in 1KG space",
    color="Super population")
# print(p)
    # ggsave(paste0("ukbb_ref.PC",i,"_PC",i+1,"_by_superpop.pdf"), p, width=7.5, height=6)
    ggsave(paste0("ukbb_ref.PC",i,"_PC",i+1,"_by_superpop.png"), p, width=7.5, height=6)
}



# Use RF to predict ancestry based on top 6 PCs
# RF function to predict ancestry using PCs:
pop_forest <- function(training_data, data, ntree=100, seed=42, pcs=1:6) {
  set.seed(seed)
  form <- formula(paste('as.factor(known_pop) ~', paste0('PC', pcs, collapse = ' + ')))
  forest <- randomForest(form,
                        data = training_data,
                        importance = T,
                        ntree = ntree)
  print(forest)
  fit_data <- data.frame(predict(forest, data, type='prob'), sample = data$sample)
  fit_data %>%
    gather(predicted_pop, probability, -sample) %>%
    group_by(sample) %>%
    slice(which.max(probability))
}

# Prep data
# Explanation:
# `trdat` and `tedat` are training and testing data. training data has to have a 
# column `known_pop` and `PC1` to `PC10` or so. `tedat` is expected to have a 
# column `sample` which is just the sample ID, and also PC columns.


trdat <- pca_comb %>%
  filter(superpop != 'UKBB') %>%
  select(superpop, PC1:PC20) %>%
  rename(known_pop = superpop)
trdat$known_pop <- as.character(trdat$known_pop)

tedat <- pca_comb %>%
  filter(superpop=='UKBB') %>%
  select(IID, PC1:PC20)
names(tedat)[1] <- "sample"
tedat$sample <- as.character(tedat$sample)


#Prediction results:
# Make prediction
pop_pred <- as.data.frame(pop_forest(training_data = trdat, data = tedat))
# when number of PCs (npc = 6)
# No. of variables tried at each split: 2

# Ntree=100
#         OOB estimate of  error rate: 0.36%
# Confusion matrix:
#     AFR AMR EAS EUR SAS class.error
# AFR 657   4   0   0   0 0.006051437
# AMR   3 342   0   2   0 0.014409222
# EAS   0   0 504   0   0 0.000000000
# EUR   0   0   0 503   0 0.000000000
# SAS   0   0   0   0 489 0.000000000

# Ntree=500 (same as 1000)
#         OOB estimate of  error rate: 0.24%
# Confusion matrix:
#     AFR AMR EAS EUR SAS class.error
# AFR 658   3   0   0   0 0.004538578
# AMR   2 344   0   1   0 0.008645533
# EAS   0   0 504   0   0 0.000000000
# EUR   0   0   0 503   0 0.000000000
# SAS   0   0   0   0 489 0.000000000

# Overall ancestry assignment
table(pop_pred$predicted_pop)
#  AFR    AMR    EAS    EUR    SAS
# 9586   4994   2563 460919  10315

# summary(pop_pred$probability)
# dim(pop_pred %>% filter(probability<0.5))
# dim(pop_pred %>% filter(probability<0.9))

summary(pop_pred$probability[pop_pred$predicted_pop=="EUR"])
 #   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
 # 0.2600  1.0000  1.0000  0.9939  1.0000  1.0000
summary(pop_pred$probability[pop_pred$predicted_pop=="AFR"])
 #   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
 # 0.2700  1.0000  1.0000  0.9542  1.0000  1.0000
summary(pop_pred$probability[pop_pred$predicted_pop=="AMR"])
 #    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
 # 0.2600  0.4200  0.5100  0.5628  0.6100  1.0000
summary(pop_pred$probability[pop_pred$predicted_pop=="EAS"])
 #   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
 # 0.2600  0.9800  1.0000  0.9456  1.0000  1.0000
summary(pop_pred$probability[pop_pred$predicted_pop=="SAS"])
 #   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
 # 0.2900  0.9900  1.0000  0.9523  1.0000  1.0000


sum(pop_pred$probability[pop_pred$predicted_pop=="AMR"] > 0.51) #2348

pop_pred.sub <- pop_pred %>%
  filter(probability >= 0.9)
table(pop_pred.sub$predicted_pop)
# when pred_prob >= 0.9 (N=473505)
#  AFR   AMR   EAS    EUR    SAS
# 8206  705   2278 453344   9290

# when pred_prob > 0.8 (N=476332)
#  AFR   AMR   EAS    EUR    SAS
# 8647   811  2293 455207   9374

# when pred_prob > 0.5 (N=484710)
#  AFR   AMR   EAS    EUR    SAS
# 9446  2578  2498 460320   9868

nrow(pop_pred) - nrow(pop_pred.sub)
# -> 14554 individuals's ancestry unclassified


##########
# Try also: nearest neighbor KNN
library(class)
knn_res <- knn(trdat[,c(2:7)], tedat[,c(2:7)], trdat$known_pop, k=20, prob=TRUE) #output is a vector of predicted pouplation labels, based on top 6 PCs
table(knn_res)
  #  AFR    AMR    EAS    EUR    SAS
  # 9332   1335   2747 464504  10459
table(knn_res[attributes(knn_res)$prob > 0.9])
  #  AFR    AMR    EAS    EUR    SAS
  # 8690    976   2676 462917  10335
##########





# Plot pca of UKBB, based on predicted super population (prob>0.9)
pca_ukb$IID <- as.character(pca_ukb$IID)
# pop_pred.sub$sample <- as.character(pop_pred.sub$sample)
# pca_ukb.pred <- merge(pca_ukb, pop_pred.sub, by.x="IID", by.y="sample")

pop_pred$predicted_pop[pop_pred$probability < 0.9] <- "Other"
pop_pred$sample <- as.character(pop_pred$sample)
pca_ukb.pred <- merge(pca_ukb, pop_pred, by.x="IID", by.y="sample")
pca_ukb.pred$predicted_pop <- factor(pca_ukb.pred$predicted_pop, levels=c("Other","AFR","AMR","EAS","EUR","SAS"))

# PC plots: ukbb samples only, colored by predicted super population/ancestry (prob > 0.9)
for(i in 1:10){
  p = ggplot(pca_ukb.pred, aes_string(x=paste0('PC',i), y=paste0('PC',i+1))) +
    geom_point(aes(color=predicted_pop), alpha=0.65, size=0.97) +
    scale_color_manual(values = c("grey25",pal_d3("category20")(20)[1:5])) +
    theme_bw() +
    guides(color = guide_legend(override.aes = list(size=2))) +
    labs(x=paste0('PC',i), y=paste0('PC',i+1), title="UKBB samples (in 1KG space)",
         color="Predicted population")
  # print(p)
  ggsave(paste0("ukbb_ref.PC",i,"_PC",i+1,"_by_predPop0.9.png"),
         p, width=7.75, height=6)
}



# Export predicted ancestral groups for study samples
# Save predicted population labels and EUR family and individual IDs
write.table(pca_ukb.pred, "ukbb_ref.projected.PC.w.pop.pred.tsv", quote=F,
  col.names=T, row.names=F, sep="\t")

# write.table(pca_ukb.pred, "ukbb_ref.projected.PC.w.pop.pred-2.tsv", quote=F,
#   col.names=T, row.names=F, sep="\t")




##### Visualize UKBB projected PCs with pre-calculated UKBB PCs: /data/tge/Tian/UKBB_full/genomics/qc/ukb_sqc_v2.txt, columns 26-45
eid <- fread("/data/tge/Tian/UKBB_full/genomics/cal_raw/ukb_cal_chr22_v2.fam", h=F)
pca_ukb.prec <- fread("/data/tge/Tian/UKBB_full/genomics/qc/ukb_sqc_v2.txt", h=F)
pca_ukb.prec <- pca_ukb.prec[,c(26:45)]
pca_ukb.prec <- cbind(eid$V1, pca_ukb.prec)
names(pca_ukb.prec) <- c("eid", paste0("PC",1:20))

pca_ukb.pred <- fread("ukbb_ref.projected.PC.w.pop.pred.tsv", h=T)
head(pca_ukb.pred,3)

pca_ukb.prec <- pca_ukb.prec[order(pca_ukb.prec$eid),]
pca_ukb.pred <- pca_ukb.pred[order(pca_ukb.pred$IID),]
sum(pca_ukb.prec$eid != pca_ukb.pred$IID)

pca_ukb.prec <- as.data.frame(pca_ukb.prec)
pca_ukb.pred <- as.data.frame(pca_ukb.pred)

for(i in 1:20){
  cat("PC",i,"\n")
  print(cor(pca_ukb.prec[,paste0("PC",i)], pca_ukb.pred[,paste0("PC",i)]))
}
# PC 1
# [1] 0.9154482
# PC 2
# [1] 0.5153534
# PC 3
# [1] -0.6791411
# PC 4
# [1] 0.174966
# PC 5
# [1] -0.06879836
# PC 6
# [1] -0.04764079


pca_ukb.pred$predicted_pop <- factor(pca_ukb.pred$predicted_pop, levels=c("Other","AFR","AMR","EAS","EUR","SAS"))

pca_ukb.prec$predicted_pop <- pca_ukb.pred$predicted_pop

# PC plots: ukbb samples only, colored by predicted super population/ancestry (prob > 0.9)
for(i in 1:10){
  p = ggplot(pca_ukb.prec, aes_string(x=paste0('PC',i), y=paste0('PC',i+1))) +
    geom_point(aes(color=predicted_pop), alpha=0.65, size=0.97) +
    scale_color_manual(values = c("grey25",pal_d3("category20")(20)[1:5])) +
    theme_bw() +
    guides(color = guide_legend(override.aes = list(size=2))) +
    labs(x=paste0('PC',i), y=paste0('PC',i+1), title="UKBB precalculated PCs",
         color="1KG-predicted population")
  # print(p)
  ggsave(paste0("ukbb.precalPC",i,"_PC",i+1,"_by_predPop0.9.png"),
         p, width=7.75, height=6)
}
