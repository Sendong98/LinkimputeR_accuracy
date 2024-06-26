#!/bin/bash
#SBATCH --job-name=linkimputer
#SBATCH --partition=regularlo
#SBATCH --nodes=1
#SBATCH --time=80:00:00
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=10G
#SBATCH --output=Linkimputer-%j.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user=s.dong@rug.nl

module load BCFtools

#select high coverage samples
bcftools view -S downsample_id_true.txt chromosome1_selected_site_biallelic.vcf.gz -Oz -o chromosome1_selected_site_biallelic_22_high.vcf.gz
bcftools index chromosome1_selected_site_biallelic_22_high.vcf.gz

#the rest samples
bcftools view -S ^downsample_id_true.txt chromosome1_selected_site_biallelic.vcf.gz -Oz -o chromosome1_selected_site_biallelic_other.vcf.gz
bcftools index chromosome1_selected_site_biallelic_other.vcf.gz

#add 10% missing values to high coverage samples
bcftools +setGT chromosome1_selected_site_biallelic_22_high.vcf.gz -Oz -o chromosome1_selected_site_biallelic_22_high_missing10.vcf.gz -- -t r:0.1 -n .
bcftools index chromosome1_selected_site_biallelic_22_high_missing10.vcf.gz

#merge missing part and rest
bcftools merge chromosome1_selected_site_biallelic_22_high_missing10.vcf.gz chromosome1_selected_site_biallelic_other.vcf.gz -Oz -o chromosome1_selected_site_biallelic_with_22_missing10.vcf.gz
bcftools index chromosome1_selected_site_biallelic_with_22_missing10.vcf.gz

#using linkimputer to perform imputation
java -Xmx200G -jar LinkImputeR.jar -s accuracy.ini

#sometimes, the output format is not correct and needs to be converted
#bcftools view chromosome1_selected_site_biallelic_with_22_missing10_imputed.vcf.gz -Oz -o chromosome1_selected_site_biallelic_with_22_missing10_imputed_new.vcf.gz

#subset non-missing sites in the original top 22 individual SNPs (chromosome1_selected_site_biallelic_22_high.vcf.gz)
bcftools view -g ^miss chromosome1_selected_site_biallelic_22_high.vcf.gz -Oz -o chromosome1_selected_site_biallelic_22_high_non_missing.vcf.gz
bcftools index chromosome1_selected_site_biallelic_22_high_non_missing.vcf.gz

#get the sites list to compare
bcftools query -f '%CHROM\t%POS\n' chromosome1_selected_site_biallelic_22_high_non_missing.vcf.gz > non_missing_sites.txt

#subset high coverage samples and non-missing sites in the original dataset to compare
bcftools view -S downsample_id_true.txt -T non_missing_sites.txt chromosome1_selected_site_biallelic_with_22_missing10_imputed.vcf.gz -Oz -o chromosome1_selected_site_biallelic_with_22_missing10_imputed_high.vcf.gz
bcftools index chromosome1_selected_site_biallelic_with_22_missing10_imputed_high.vcf.gz

#get bcftools gtcheck compare file
awk '{print $1,$1}' downsample_id_true.txt > gtcheck22.txt

#check discordance using bcftools gtcheck
bcftools gtcheck -P gtcheck22.txt -g chromosome1_selected_site_biallelic_with_22_missing10_imputed_high.vcf.gz chromosome1_selected_site_biallelic_22_high_non_missing.vcf.gz -e 0 > discordance_chr1.txt

