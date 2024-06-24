#!/bin/bash
#SBATCH --job-name=linkimputer
#SBATCH --partition=regularlo
#SBATCH --nodes=1
#SBATCH --time=80:00:00
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=30G
#SBATCH --output=Impute-%j.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user=s.dong@rug.nl

#select high coverage samples
bcftools view -S ../../LinkImputeR/downsample_id_true.txt chromsome1_selected_site_biallellic.vcf.gz -Oz -o chromsome1_selected_site_biallellic_22_high.vcf.gz

#the rest samples
bcftools view -S ^../../LinkImputeR/downsample_id_true.txt chromsome1_selected_site_biallellic.vcf.gz -Oz -o chromsome1_selected_site_biallellic_other.vcf.gz

#add 10% missing values to high coverage samples
bcftools +setGT chromsome1_selected_site_biallellic.vcf.gz -Oz -- -t r:0.1 -n . -o test.vcf.gz

#merge missing part and rest
bcftools merge chromsome1_selected_site_biallellic_22_high_missing10.vcf.gz chromsome1_selected_site_biallellic_other.vcf.gz -Oz -o chromsome1_selected_site_biallellic_with_22_missing10.vcf.gz

#using linkimputer to perform imputation
java -Xmx200G -jar LinkImputeR.jar -s accuracy.ini

#check accuracy
