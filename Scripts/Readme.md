#Methods
##Method 1: imputeqc
imputedqc package provides masking and accuracy checking function
Mask 10% proportion of missing genotypes, only one file to be generated
`Rscript mask_imputeqc.R -n 1 -p 0.1 -o chr1_missing10 chromosome1_selected_site_biallelic.vcf.gz
bcftools view chr1_missing10.m1.vcf -Oz -o chr1_missing10.m1.vcf.gz`
Check the discordance
`Rscript Imputeqc.R`

##Method 2: linkimputer.sh
