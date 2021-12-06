#!/bin/bash

#SBATCH -A m2_jgu-funcpoll
#SBATCH -J muscle

#SBATCH -p devel
#SBATCH -N 1
#SBATCH -n 40
#SBATCH --cpu-per-task 1
#SBATCH -t 04:00:00

#SBATCH -v
#SBATCH -o muscle.%j.out
#SBATCH -e muscle.%j.err

module purge
module load bio/MUSCLE/3.8.1551-GCC-9.3.0

PROJECT_DIR=/lustre/project/m2_jgu-funcpoll
IN_DIR=${PROJECT_DIR}/phylogenetic_analysis/bee_proteomes/results/orthofinder/Results_Oct19/Single_Copy_Orthologue_Sequences
OUT_DIR=${PROJECT_DIR}/phylogenetic_analysis/fasta/04_sco_protein_alignments_v2

mkdir -p ${OUT_DIR}

for ORTHOGROUP_FASTA in $(ls ${IN_DIR});
do
	echo ${ORTHOGROUP_FASTA}
	srun -n 1 -c 1 muscle -in ${IN_DIR}/${ORTHOGROUP_FASTA} -out ${OUT_DIR}/${ORTHOGROUP_FASTA} &
done

wait
