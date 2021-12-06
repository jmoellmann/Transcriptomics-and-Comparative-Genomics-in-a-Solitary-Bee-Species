#!/bin/bash

#SBATCH -A m2_jgu-funcpoll
#SBATCH -J star_after_trimming

#SBATCH -p devel
#SBATCH -n 27
#SBATCH -c 5
#SBATCH -N 4
#SBATCH --mem-per-cpu 2G
#SBATCH -t 04:00:00

#SBATCH -v
#SBATCH -o star_after_trimming.%j.out
#SBATCH -e star_after_trimming.%j.err

module purge
module load bio/STAR/2.7.3a-GCC-9.3.0

PROJECT_DIR=/lustre/project/m2_jgu-funcpoll
EXPERIMENT_DIR=${PROJECT_DIR}/data/datasets/PRJEB9366
INDEX_DIR=${PROJECT_DIR}/data/reference/STAR_index_Bter_1.0

FASTQ_DATA_DIR=${EXPERIMENT_DIR}/fastp-trimmed
STAR_OUT_DIR=${PROJECT_DIR}/results/PRJEB9366/after_trimming/alignment

mkdir -p ${STAR_OUT_DIR}

for FASTQ_FILE in $(ls ${FASTQ_DATA_DIR});
do
	SAMPLE=${FASTQ_FILE%"_trimmed.fastq.gz"}

	if [ ! -f  ${STAR_OUT_DIR}/${SAMPLE}_Log.final.out ]
	then

		IN_FILE=${FASTQ_DATA_DIR}/${FASTQ_FILE}

		srun -n 1 -c 5 --exclusive \
		STAR --runThreadN 5 --genomeDir ${INDEX_DIR} \
		--readFilesIn ${IN_FILE} --readFilesCommand gunzip -c \
		--genomeLoad NoSharedMemory --outFileNamePrefix ${STAR_OUT_DIR}/${SAMPLE}_ \
		--outSAMtype BAM Unsorted \
		--limitBAMsortRAM 8500000000 &

	fi

done
wait
