#!/bin/bash

#SBATCH --mail-type=end,fail
#SBATCH --job-name="MAKER"
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --time=96:00:00
#SBATCH --mem=25G
#SBATCH --partition=pall


NAME="10klist.fa"

module add SequenceAnalysis/GenePrediction/augustus/3.3.3.1
cp -r /software/SequenceAnalysis/GenePrediction/augustus/3.3.3.1/config ./augustus_config
export AUGUSTUS_CONFIG_PATH="./augustus_config/"

module add SequenceAnalysis/GenePrediction/maker/2.31.9

maker -fix_nucleotides maker_opts.ctl maker_bopts.ctl maker_exe.ctl

#then merge gffs
gff3_merge -n -d ${NAME}.maker.output/${NAME}_master_datastore_index.log

#results: 10klist.all.gff

#get protein fasta
find ${NAME}.maker.output/*/*/* -name *maker.proteins.fasta -exec cat {} >> allprot.fa \;

