#!/usr/bin/env bash
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1G
#SBATCH --time=3:00:00
#SBATCH --job-name=yeast_SNPcalling
#SBATCH --mail-user=aimee.freiberg@students.unibe.ch
#SBATCH --mail-type=begin,end,fail

####################
###	SNP calling ###
###################

#generate variables
USER=afreiberg

#load modules
source /data/users/lfalquet/SBC07107_21/scripts/module.sh

#create and go to the TP directory
mkdir /data/users/${USER}/Bioinformatics/snp_calling
cd /data/users/${USER}/Bioinformatics/snp_calling

#link the reference genome and the reads locally
ln -s /data/users/lfalquet/SBC07107_21/reference/R64-1-1.104.fa R64-1-1.104.fa
ln -s /data/users/lfalquet/SBC07107_21/reference/R64-1-1.104.gtf R64-1-1.104.gtf
ln -s /data/users/${USER}/Bioinformatics/qc/*.bam .

while read line; do

# call SNPs for each strain (each bam file) and index the vcf.gz  -- D10_S108
bcftools mpileup -Ou --threads 8 -f R64-1-1.104.fa ${line}.bam | bcftools call -vc -Oz --threads 8 -o ${line}.vcf.gz
tabix ${line}.vcf.gz

done < strains_name.txt

#or much faster merged with bcftools once you have all vcf files and their indexes
bcftools merge -O v `ls *.vcf.gz` > D10_15all.vcf

#annotate vcf file using snpEff
#first download and install snpEff locally
wget https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip
unzip snpEff_latest_core.zip

# get the database for the yeast genome (warning in snpEff 99=104).
java -Xmx4g -jar ./snpEff/snpEff.jar download R64-1-1.99

#annotate the VCF file
java -Xmx20g -jar ./snpEff/snpEff.jar -no-upstream -no-downstream R64-1-1.99 D10_15all.vcf > D10_15all_annot.vcf

#remove synonymous and intergenic variants
cat D10_15all_annot.vcf | java -jar ./snpEff/SnpSift.jar filter "(( ANN[*].EFFECT != 'synonymous_variant') & ( ANN[*].EFFECT != 'intergenic_region'))" > D10_15all_coding.vcf

# keep only variants that are found in less than 4 strains
cat D10_15all_coding.vcf | java -jar ./snpEff/SnpSift.jar filter "((countVariant() < 4))" > D10_15all_coding_max3.vcf

# Tip: modify before importing into excel
perl -pe 's/;ANN=/;\tANN=/g' D10_15all_coding_max3.vcf > D10_15all_coding_max3_forexcel.vcf
perl -i -pe 's/FORMAT\t/FORMAT\t\t/' D10_15all_coding_max3_forexcel.vcf