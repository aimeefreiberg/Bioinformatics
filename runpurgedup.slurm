#!/bin/bash

#SBATCH --mail-type=fail
#SBATCH --job-name="hifiasm"
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --time=24:00:00
#SBATCH --mem=400G
#SBATCH --partition=pcoursea

NAME=$1 #fasta/q file of reads

#no purging (avoid haploids)
#/home/lfalquet/softs/hifiasm -o asm_hifi_tsetse -t40 -l0 ${NAME} 2> hifiasm.log
#perl -ane 'if ($F[0]eq"S") { print ">".$F[1]."\n".$F[2]."\n"}' asm_hifi_tsetse.bp.p_ctg.gfa > asm_hifi_tsetse.p_ctg.fa

cp /data/users/lfalquet/SBC07107_21/scripts/pb.fofn .
/home/lfalquet/softs/purge_dups/scripts/pd_config.py -n configpurge.json ${NAME} pb.fofn 
perl -i -pe 's/(normal|long)/pcoursea/g' configpurge.json

#ignore error
#edit config.json to accommodate slurm or copy from scripts
#cp /data/users/lfalquet/SBC07107_21/scripts/pb.fofn .
#cp /data/users/lfalquet/SBC07107_21/scripts/config.json .
module add UHTS/Analysis/minimap2/2.17
python3 /home/lfalquet/softs/purge_dups/scripts/run_purge_dups.py -p slurm configpurge.json /home/lfalquet/softs/purge_dups/bin Tsetse

