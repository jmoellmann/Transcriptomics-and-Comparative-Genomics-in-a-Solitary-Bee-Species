#!/bin/bash

#SBATCH -A m2_jgu-funcpoll
#SBATCH -J get_cds_for_msa

#SBATCH -p devel
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -t 04:00:00

#SBATCH -v
#SBATCH -o get_cds_for_msa.%j.out
#SBATCH -e get_cds_for_msa.%j.err

module purge
module load bio/HTSeq/0.13.5-foss-2020a-Python-3.8.2

srun -n 1 -N 1 python /lustre/project/m2_jgu-funcpoll/scripts/PRJNA285788/phylogenetic_analysis/general/get_cds_for_msa.py &

wait

