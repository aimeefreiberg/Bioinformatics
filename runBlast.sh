#!/usr/bin/env bash
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1G
#SBATCH --time=3:00:00
#SBATCH --job-name=yeast_SNPcalling
#SBATCH --mail-user=aimee.freiberg@students.unibe.ch
#SBATCH --mail-type=begin,end,fail

######################################################
### Blast query for hematophagy genes and proteins ###
######################################################

module add Blast/ncbi-blast/2.10.1+
makeblastdb -in allprot.fa -dbtype prot
blastp -query proteins_aegypti.fasta -db allprot.fa
