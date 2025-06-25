######################### Input data #########################

###Atopic Dermatitis:
#Study, GCST90244788:
atopic_dermatitis.gwas.gcst90244788.tsv 

#Study, GSE121212, paired:
atopic_dermatitis.rnaseq.gse121212.tsv

#Fisher integrated:
atopic_dermatitis.integrated.tsv



### Focal Epilepsy:
#Study, Kjaer_Guelfi:
focal_epilepsy.rnaseq.kjaer_guelfi_consensus.txt



### Ulcerative Colitis:
#Study, GSE66407, GSE109142 and GSE166925:
ulcerative_colitis.rnaseq.integrated.tsv


#Study, GWAS Finngen and GWAS 2015:
ulcerative_colitis.gwas.integrated.tsv


#Fisher integrated:
ulcerative_colitis.integrated.tsv


######################### Gold standard #########################

# gold standard positive negative balanced
atopic_dermatitis.gold_standard.balanced.tsv  
focal_epilepsy.gold_standard.balanced.tsv  
ulcerative_colitis.gold_standard.balanced.tsv

# gold standard not positive negative balanced
atopic_dermatitis.gold_standard.tsv           
focal_epilepsy.gold_standard.tsv           
ulcerative_colitis.gold_standard.tsv



######################### Protein embedding #########################

9606.node2vec64.h5

