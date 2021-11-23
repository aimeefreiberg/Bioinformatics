#!/bin/bash

#SBATCH --mail-type=fail
#SBATCH --job-name="asm"
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --time=24:00:00
#SBATCH --mem=40G
#SBATCH --partition=pcoursea

source /data/users/lfalquet/SBC07107_21/scripts/module.sh

XX=$1

#create and go to the TP directory
mkdir asm_${XX}
cd asm_${XX}

#link the reference genome and the reads locally
ln -s /data/users/lfalquet/SBC07107_21/hifireads/${XX}.fastq . 
#verify that your files are in your directory (ls -l)

#check quality of your data with fastqc (don't do that as it requires too much RAM)
#fastqc -t 4 ${XX}.fastq

#optional view results from your public-html or download the html file to your local machine
#scp student@168.176.54.15:?? ~/Desktop/
#(use winscp on windows)

#clean reads
fastp -i ${XX}.fastq -W 10 -5 -3 -M 10 -l 1000 -G -A -w 4 -o ${XX}_clean.fastq 

#assemble with hifiasm
hifiasm -o asm_hifiasm -t40 ${XX}_clean.fastq 2> hifiasm.log
perl -ane 'if ($F[0]eq"S") { print ">".$F[1]."\n".$F[2]."\n"}' asm_hifiasm.p_ctg.gfa > asm_hifiasm.p_ctg.fa



