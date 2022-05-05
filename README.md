# ukbb-pca

This directory details the steps performed to infer genetic ancestry of the UKBB samples (used in the PRS-CSx project).


## PRS-CSx: population assignment

### Methods

* Intersect 1KG with UKBB genotyped SNPs
* Based on the overlapping set of SNPs, 
	* Remove strand ambiguous SNPs
	* Remove long-range LD regions
	* Keep those with call rate > 0.98 (in 1KG)
	* Keep those with MAF > 0.05 (in 1KG)
	* Keep autosomal SNPs only (in 1KG)
* Perform LD pruning (plink --indep-pairwise 100 50 0.2) to obtain a list of independent high-quality markers
	* 149501 variants left
* Run PCA in 1KG
* Project PC loadings onto UKBB samples
	* With scale appropriately adjusted (more details on this)
* Classify ancestry for UKBB samples
	* Random Forest classifier based on top 6 PCs with 1KG as the training set
	* Assign population label based on a prediction probability > 0.9
* For each predicted population, perform sample QC
	* Remove the following samples:
		* Sex check (submitted != inferred)
		* Missingness/Heterozygosity outliers
		* Sex chromosome aneuploidy
		* Related individuals (not in “used.in.pca.calculation”)
	* Post-QC:
		* For EUR, separate those in Neale lab GWAS (training) vs. not (testing)
		* For the other ancestral groups, save a copy for each (training)
* Run PCA for each population in UKBB
	* erform population-specific LD pruning based on a list of high-quality markers (same marker QC as described above)
	* Run PCA


### Summary

| Super population    |  Predicted N   |  Overlap with Neale Lab GWAS   |  Predicted N, post-QC*   |  In Neale Lab: YES   |  In Neale Lab: NO  |
| --- | -----: | -----: | -----: | -----: | -----: |   
| AFR | 8,206 | 0 | 7,507 | 0 | 7,507 |
| AMR | 705 | 0 | 687 | 0 | 687 |
| EAS | 2,278 | 0 | 2,181 | 0 | 2,181 |
| EUR | 453,344 | 361,035 | 375,120 | 361,035 | 14,085 |
| SAS | 9,290 | 0 | 8,412 | 0 | 8,412 |
| Other | 14,554 | 159 | 12,905 | 159 | 12,746 |
| Total | 488,377 | 361,194 | 406,812 | 361,194 | 45,618 |

