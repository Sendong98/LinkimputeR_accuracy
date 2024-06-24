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

#select high coverage samples
bcftools view -S downsample_id_true.txt chromsome1_selected_site_biallellic.vcf.gz -Oz -o chromsome1_selected_site_biallellic_22_high.vcf.gz

#the rest samples
bcftools view -S ^downsample_id_true.txt chromsome1_selected_site_biallellic.vcf.gz -Oz -o chromsome1_selected_site_biallellic_other.vcf.gz

#add 10% missing values to high coverage samples
bcftools +setGT chromsome1_selected_site_biallellic_22_high.vcf.gz -Oz -- -t r:0.1 -n . -o chromsome1_selected_site_biallellic_22_high_missing10.vcf.gz

#merge missing part and rest
bcftools merge chromsome1_selected_site_biallellic_22_high_missing10.vcf.gz chromsome1_selected_site_biallellic_other.vcf.gz -Oz -o chromsome1_selected_site_biallellic_with_22_missing10.vcf.gz

#using linkimputer to perform imputation
java -Xmx200G -jar LinkImputeR.jar -s accuracy.ini

#Subset high coverage samples to compare
bcftools view -S downsample_id_true.txt chromsome1_selected_site_biallellic_with_22_missing10_imputed.vcf.gz -Oz -o chromsome1_selected_site_biallellic_with_22_missing10_imputed_high.vcf.gz

#check accuracy

