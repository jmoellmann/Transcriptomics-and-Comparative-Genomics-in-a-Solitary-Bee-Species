#!/bin/bash

#SBATCH -A m2_jgu-funcpoll
#SBATCH -J hyphy

#SBATCH -p parallel
#SBATCH -N 6
#SBATCH -n 240
#SBATCH -c 1
#SBATCH --mem-per-cpu=1G
#SBATCH -t 04:00:00

#SBATCH -v
#SBATCH -o hyphy.%j.out
#SBATCH -e hyphy.%j.err

module purge

PROJECT_DIR=/lustre/project/m2_jgu-funcpoll
HYPHY_DIR=${PROJECT_DIR}/phylogenetic_analysis/src/hyphy/

SINGLE_COPY_ORTHOLOGUES=${PROJECT_DIR}/phylogenetic_analysis/fasta/Orthogroups_SingleCopyOrthologues.txt
TREES_DIR=${PROJECT_DIR}/phylogenetic_analysis/bee_proteomes/results/orthofinder/Results_Oct19/Gene_Trees
ALIGNMENTS_DIR=${PROJECT_DIR}/phylogenetic_analysis/fasta/06_cds_alignments_v2
OUT_DIR=${PROJECT_DIR}/phylogenetic_analysis/fasta/hyphy_results

mkdir -p ${OUT_DIR}
cd ${HYPHY_DIR}

mkdir ${OUT_DIR}/log

for ORTHOGROUP in $(cat ${SINGLE_COPY_ORTHOLOGUES});
do
	touch ${OUT_DIR}/log/${ORTHOGROUP}_log.txt
	srun -n 1 -c 1 --exclusive ./hyphy absrel --alignment ${ALIGNMENTS_DIR}/${ORTHOGROUP}.fa --tree ${TREES_DIR}/${ORTHOGROUP}_tree.txt \
	--output ${OUT_DIR}/${ORTHOGROUP}.txt > ${OUT_DIR}/log/${ORTHOGROUP}_log.txt &

done

wait
