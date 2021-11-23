#!/bin/bash

#SBATCH --mail-type=end,fail
#SBATCH --job-name="bwa"
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --time=3:00:00
#SBATCH --mem=25G
#SBATCH --partition=pcoursea

XX=$1

#create and go to the TP directory
mkdir bwa$1
cd bwa$1

source /data/users/lfalquet/SBC07107_21/scripts/module.sh

#link the reference genome and the reads locally
ln -s /data/users/lfalquet/SBC07107_21/reference/R64-1-1.104.fa .
ln -s /data/users/lfalquet/SBC07107_21/rawdata/${XX}_R1*.fastq.gz ${XX}_1.fastq.gz
ln -s /data/users/lfalquet/SBC07107_21/rawdata/${XX}_R2*.fastq.gz ${XX}_2.fastq.gz

#verify that your files are in your directory (ls -l)

fastqc -t 2 ${XX}_1.fastq.gz ${XX}_2.fastq.gz

#with trimmotatic
#trimmomatic PE -threads 4 -phred33 ${XX}_1.fastq.gz ${XX}_2.fastq.gz ${XX}_1trim.fastq.gz ${XX}_1unpaired.fastq.gz ${XX}_2trim.fastq.gz ${XX}_2unpaired.fastq.gz LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:130

#with fastp
fastp -i ${XX}_1.fastq.gz -I ${XX}_2.fastq.gz -o ${XX}_1trim.fastq.gz -O ${XX}_2trim.fastq.gz -j ${XX}_fastp.json -h ${XX}_fastp.html --thread 4 --trim_poly_g -l 150;


#index reference for bwa
bwa index R64-1-1.104.fa

#map reads onto reference
bwa mem -t 4 -M R64-1-1.104.fa ${XX}_1trim.fastq.gz ${XX}_2trim.fastq.gz > ${XX}.sam

#index reference for samtools
samtools faidx R64-1-1.104.fa

#convert sam to bam
samtools view -b -@ 4 -t R64-1-1.104.fa.fai ${XX}.sam > ${XX}_unsorted.bam

#sort bam
samtools sort -@ 4 -o ${XX}.bam ${XX}_unsorted.bam
rm ${XX}_unsorted.bam  
rm ${XX}.sam

#index bam
samtools index ${XX}.bam

