# LinkimputeR_accuracy
These script are using for the accuracy test of LinkimputeR. 
Genomic data: Seychelles warbler chromosome 1
Accuracy check method: downsampling the 22 top coverage individuals (coverage over 8*) and impute with other individuals (2006 in total), compare with true dataset and imputed dataset by using imputeqc and bcftools gtcheck
