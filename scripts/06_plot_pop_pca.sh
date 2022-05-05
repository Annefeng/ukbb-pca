
################################
# PRS-CSx: population assignment
################################

module load R/3.4.0
source ~/R-3.4.0-ownlib.bash
R

setwd("/data/js95/yfeng/projects/prscsx_pca/results")

library(data.table)
library(randomForest)
library(tidyverse)
library(ggsci)
library(RColorBrewer)

pop="eur.non.neale2"
for( pop in c("eur.non.neale","afr","amr","eas","sas") ){

    print(pop)
    # Read in data
    pca <- fread(paste0("ukbb_sqc_",pop,"_predProb0.9.autsnp.geno02.maf05.pca.eigenvec"), h=T)
    pca$pop <- substring(pop,1,3)

    # Plot pca of UKBB predicted population (prob>0.9)
    for(i in 1:5){
      p = ggplot(pca, aes_string(x=paste0('PC',i), y=paste0('PC',i+1))) +
        geom_point(aes(color=pop), alpha=0.65, size=0.97) +
        scale_color_manual(values = pal_d3("category20")(20)[1]) +
        theme_bw() +
        guides(color = guide_legend(override.aes = list(size=2))) +
        labs(x=paste0('PC',i), y=paste0('PC',i+1), title=paste0("UKBB ",toupper(pop)," samples"))
      # print(p)
      ggsave(paste0("ukbb_sqc_",pop,"_predPop0.9_PC",i,"_PC",i+1,".png"),
             p, width=6.75, height=6)
  }

}




##### EUR with 1kg
pca.1kg.eur <- read.table("ukbb.1kg_eur.pca.eigenvec", header=T, sep="\t", stringsAsFactors=F)
pca.1kg.eur$pop <- pca.1kg.eur$FID
pca.1kg.eur$pop[pca.1kg.eur$pop==pca.1kg.eur$IID] <- "UKBB"
table(pca.1kg.eur$pop)

pca.1kg.eur$pop <- factor(pca.1kg.eur$pop,levels=c("UKBB","CEU","FIN","GBR","IBS","TSI"))

for(i in 1:5){
  p = ggplot(pca.1kg.eur, aes_string(x=paste0('PC',i), y=paste0('PC',i+1))) +
    geom_point(aes(color=pop), alpha=0.55) +
    scale_color_manual(values = pal_d3("category20")(20)[1:6]) + 
    theme_bw() +
    guides(color = guide_legend(override.aes = list(size=2))) +
    labs(x=paste0('PC',i), y=paste0('PC',i+1), title="UKBB + 1KG EUR samples",
         color="1KG populations")
  # print(p)
  ggsave(paste0("ukbb_sqc.1kg_eur_PC",i,"_PC",i+1,"_by_pop.png"),
       p, width=7.5, height=6)
}




###
pop="eur.non.neale"
pca <- fread(paste0("ukbb_sqc_",pop,"_predProb0.9.autsnp.geno02.maf05.pca.eigenvec"), h=T)
pca$pop <- substring(pop,1,3)

sr_race <- fread("/data/js95/yfeng/projects/ukbb_aoo/misc_data/ukb32568_20190220_eid_21000.tsv", h=T)
names(sr_race)[2] <- "ethnicity"
table(sr_race$ethnicity)

pca_race <- merge(pca, sr_race, by.x="IID", by.y="eid")
pca$IID[!pca$IID %in% pca_race$IID] #1542021 3682334 4084900
table(pca_race$ethnicity)
  # -3   -1    1    2    4    5    6 1001 1002 1003 2001 2002 2003 2004 3004 4002
  # 57    7   38    6    2    1  357 7459   96 5859    5    3   11  151    1    1


pca_race$ethnicity <- factor(pca_race$ethnicity, levels=c("1001","1002","1003",
    "2001","2002","2003","2004","3004","4002","1","2","4","5","6","-1","-3"))
table(pca_race$ethnicity)
levels(pca_race$ethnicity) <- c("British","Irish","Any other White background",
    "White and Black Caribbean","White and Black African","White and Asian",
    "Any other mixed background","Any other Asian background","African",
    "White","Mixed","Black or Black British","Chinese","Other ethnic group",
    "Do not know","Prefer not to answer")
for(i in 1:5){
  p = ggplot(pca_race, aes_string(x=paste0('PC',i), y=paste0('PC',i+1))) +
    geom_point(aes(color=ethnicity), alpha=0.7) +
    scale_color_manual(values = pal_d3("category20")(20)) + 
    theme_bw() +
    guides(color = guide_legend(override.aes = list(size=2))) +
    labs(x=paste0('PC',i), y=paste0('PC',i+1), title="UKBB EUR samples",
         color="Reported ethnicity")
  # print(p)
  ggsave(paste0("ukbb_sqc_",pop,"_predPop0.9_PC",i,"_PC",i+1,"-by-ethnicity.png"),
       p, width=8.5, height=6)
}
