
# Intersect classified EUR with Neale lab GWAS samples

# Population lables (classified using random forest): all 488377 inds
# /data/js95/yfeng/projects/prscsx_pca/results/ukbb_ref.projected.PC.w.pop.pred.tsv

# Neale lab GWAS samples (individual IDs in Tian's app)
# /data/js95/yfeng/projects/ukbb_aoo/genetic_data/ukb32568.gwas_samples.txt



MISC=/data/js95/yfeng/projects/ukbb_aoo/misc_data

awk '$23=="EUR"{print $1"\t"$1}' /data/js95/yfeng/projects/prscsx_pca/results/ukbb_ref.projected.PC.w.pop.pred.tsv > $MISC/ukb32568_20190412_eur_genomics.fam
# -> 453344 classified EUR inds

awk 'NR==FNR{a[$1];next} !($1 in a) {print $0}' /data/js95/yfeng/projects/ukbb_aoo/genetic_data/ukb32568.gwas_samples.txt $MISC/ukb32568_20190412_eur_genomics.fam > $MISC/ukb32568_20190412_eur_genomics_nonNeale.fam
# -> 92309 EUR not in Neale
