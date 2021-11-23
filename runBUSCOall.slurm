#!/bin/bash

#SBATCH --mail-type=end,fail
#SBATCH --job-name="BUSCO"
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=3:00:00
#SBATCH --mem=30G
#SBATCH --partition=pcoursea

cd /data/users/$USER/Compare_asm

for k in `ls *.fa`
do
#y=${k%.fasta}
#sbatch /data/users/lfalquet/SBC07107_21/scripts/runBUSCO.slurm  ${y##*/}
sbatch /data/users/lfalquet/SBC07107_21/scripts/runBUSCO.slurm  ${k}
done

echo "BUSCO launched"
