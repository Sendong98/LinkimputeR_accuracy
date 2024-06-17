# LinkimputeR_accuracy
These script are using for the accuracy test of LinkimputeR. 
Genomic data: Seychelles warbler chromosome 10
Accuracy check method: downsampling the 22 top coverage individuals (coverage over 8*) and impute with other individuals (2006 in total), compare true dataset (chr10.vcf.gz) and imputed dataset (chr10_imputed.vcf.gz)

Method 1
1. Using imputeqc package: mask 10% of vcf.gz data (using make_masking_files.R), generate only 1 mask file, impute genotype with LinkimputeR, and check accuracy with Estimateqc().
