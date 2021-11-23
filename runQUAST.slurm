#!/bin/bash

#SBATCH --mail-type=end,fail
#SBATCH --job-name="QUAST"
#SBATCH --nodes=1
#SBATCH --cpus-per-task=12
#SBATCH --time=3:00:00
#SBATCH --mem=20G
#SBATCH --partition=pshort


cd /data/users/$USER

module add UHTS/Quality_control/quast/4.6.0
cp -R /data/users/lfalquet/SBC07107_21/Tsetse/Compare_asm Compare_asm

cd /data/users/$USER/Compare_asm
python /software/UHTS/Quality_control/quast/4.6.0/bin/quast.py -t 12 -o quast_out *.fa
