#!/bin/bash

#SBATCH -A m2_jgu-funcpoll
#SBATCH -J salmon_index

#SBATCH -p parallel
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --mem-per-cpu 2G
#SBATCH -t 01:00:00

#SBATCH -v
#SBATCH -o salmon_index.%j.out
#SBATCH -e salmon_index.%j.err
module purge
module load bio/Salmon/1.4.0-gompi-2020a

PROJECT_DIR=/lustre/project/m2_jgu-funcpoll
ANNOTATION_FASTA_FILE=${PROJECT_DIR}/data/reference/annotation/Apis_mellifera.Amel_HAv3.1.49.cdna.all.ncrna.plus.genome.fa
INDEX_OUT_DIR=${PROJECT_DIR}/data/reference/salmon_index
DECOYS_FILE=${PROJECT_DIR}/data/reference/annotation/decoys.txt

mkdir -p ${INDEX_OUT_DIR}

salmon index -t ${ANNOTATION_FASTA_FILE} -i ${INDEX_OUT_DIR} \
--decoys ${DECOYS_FILE}


