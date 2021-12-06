#!/bin/bash

module purge
module load bio/Pysam/0.16.0.1-GCC-9.3.0

PROJECT_DIR=/lustre/project/m2_jgu-funcpoll
BAM_DIR=${PROJECT_DIR}/results/PRJNA533756/after_trimming/quantification/kallisto
SCRIPT=${PROJECT_DIR}/scripts/PRJNA533756/annotate_bam_for_rmats.py

for SAMPLE in $(ls ${BAM_DIR});
do
	python ${SCRIPT} ${BAM_DIR}/${SAMPLE}/pseudoalignments.bam
done

