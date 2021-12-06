#!/bin/bash

module purge
module load bio/kallisto/0.46.1-foss-2019b 

PROJECT_DIR=/lustre/project/m2_jgu-funcpoll
TRANSCRIPTOME_FILE=${PROJECT_DIR}/data/reference/transcriptome/GCF_004153925.1_Obicornis_v3_rna.fna
OUT_DIR=${PROJECT_DIR}/data/reference/kallisto_index

mkdir -p ${OUT_DIR}

kallisto index ${TRANSCRIPTOME_FILE} -i ${OUT_DIR}/Osmia_bicornis_v3.kallisto_index.idx
