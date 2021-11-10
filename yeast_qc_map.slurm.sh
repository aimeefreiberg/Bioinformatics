#!/usr/bin/env bash
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1G
#SBATCH --time=3:00:00
#SBATCH --job-name=yeast_fastqualitycheck
#SBATCH --mail-user=aimee.freiberg@students.unibe.ch
#SBATCH --mail-type=begin,end,fail

########################################
###	Raw data quality check + mapping ###
########################################

#set up data
USER=afreiberg
source /data/users/lfalquet/SBC07107_21/scripts/module.sh

#link refrence genome
cd /data/users/afreiberg/Bioinformatics

ln -s /data/users/lfalquet/SBC07107_21/reference/R64-1-1.104.fa R64-1-1.104.fa
ln -s /data/users/lfalquet/SBC07107_21/reference/R64-1-1.104.gtf R64-1-1.104.gtf

mkdir /data/users/${USER}/Bioinformatics/qc

cd /data/users/afreiberg/Bioinformatics/qc

while read line; do

#copy needed file
ln -s /data/users/lfalquet/SBC07107_21/rawdata/${line}.fastq.gz ${line}.fastq.gz

#run quality check on file
fastqc -t 2 ${line}.fastq.gz
done < strains.txt

#check all quality
module add UHTS/Analysis/MultiQC/1.8
multiqc *.zip


#clean files with with fastp

while read line; do
#trim
fastp -i ${line}_R1_001.fastq.gz -I ${line}_R2_001.fastq.gz -o ${line}_1trim.fastq.gz -O ${line}_2trim.fastq.gz -j ${line}_fastp.json -h ${line}_fastp.html --thread 4 --trim_poly_g -l 150;


#quality check again
fastqc -t 2 ${line}_1trim.fastq.gz
fastqc -t 2 ${line}_2trim.fastq.gz

#index reference for bwa
bwa index R64-1-1.104.fa

#map reads onto reference - 
bwa mem -t 4 -M R64-1-1.104.fa ${line}_1trim.fastq.gz ${line}_2trim.fastq.gz > ${line}.sam

#index reference for samtools
samtools faidx R64-1-1.104.fa

#convert sam to bam
samtools view -b -@4 -tR64-1-1.104.fa.fai ${line}.sam > ${line}_unsorted.bam

#sort bam
samtools sort -@4 -o ${line}.bam ${line}_unsorted.bam
rm ${line}_unsorted.bam
rm ${line}.sam

#index bam -alignment file of all the reads to the refrence genome
samtools index ${line}.bam

done < strains_name.txt

#extract quality for all trims
module add UHTS/Analysis/MultiQC/1.8
multiqc *trim_fastqc.zip