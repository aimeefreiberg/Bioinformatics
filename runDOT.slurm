#!/bin/bash

#SBATCH --mail-type=end,fail
#SBATCH --job-name="DOT"
#SBATCH --nodes=1
#SBATCH --cpus-per-task=12
#SBATCH --time=3:00:00
#SBATCH --mem=50G
#SBATCH --partition=pcoursea

source /data/users/lfalquet/SBC07107_21/scripts/module.sh
cd /data/users/$USER/Compare_asm

nucmer -t 12 asm_hifi_tsetse.p_ctg.purged.fa asm_hifi_tsetse.p_ctg.fa -p hifiasm_purgevsraw
DotPrep.py --delta hifiasm_purgevsraw.delta --out hifiasm_purgevsraw_dot
