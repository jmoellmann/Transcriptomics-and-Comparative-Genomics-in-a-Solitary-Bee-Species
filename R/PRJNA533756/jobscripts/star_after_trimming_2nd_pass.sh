#!/bin/bash

#SBATCH -A m2_jgu-funcpoll
#SBATCH -J star_after_trimming

#SBATCH -C skylake
#SBATCH -p smp
#SBATCH -n 3
#SBATCH -c 5
#SBATCH --mem-per-cpu 2G
#SBATCH -t 12:00:00

#SBATCH -v
#SBATCH -o star_after_trimming.%j.out
#SBATCH -e star_after_trimming.%j.err

module purge
module load bio/STAR/2.7.3a-GCC-9.3.0

PROJECT_DIR=/lustre/project/m2_jgu-funcpoll
EXPERIMENT_DIR=${PROJECT_DIR}/data/datasets/PRJNA533756
INDEX_DIR=${PROJECT_DIR}/data/reference/STAR_index

for GROUP in "worker_nurse" "worker_forager" "adult_male" "virgin_queen" "mated_queen";
do
	GROUP_ACCESSION=${EXPERIMENT_DIR}/accession/${GROUP}.txt
	FASTQ_DATA_DIR=${EXPERIMENT_DIR}/fastp-trimmed/${GROUP}
	STAR_OUT_DIR=${PROJECT_DIR}/results/PRJNA533756/after_trimming/alignment/${GROUP}
	SJ_FILES=$(ls ${PROJECT_DIR}/data/reference/STAR_SJ_files/PRJNA533756/${GROUP}/*)

	mkdir -p ${STAR_OUT_DIR}

	for SAMPLE in $(cat ${GROUP_ACCESSION});
	do

		if [ ! -f  ${STAR_OUT_DIR}/${SAMPLE}_Log.final.out ]
		then

			IN_FILE=${FASTQ_DATA_DIR}/${SAMPLE}_trimmed.fastq.gz

			srun -n 1 -c 5 \
			STAR --runThreadN 5 --genomeDir ${INDEX_DIR} \
			--readFilesIn ${IN_FILE} --readFilesCommand gunzip -c \
			--genomeLoad NoSharedMemory --outFileNamePrefix ${STAR_OUT_DIR}/${SAMPLE}_ \
			--outSAMtype BAM SortedByCoordinate --quantMode GeneCounts \
			--sjdbFileChrStartEnd ${SJ_FILES} \
			--limitBAMsortRAM 8500000000 &

		fi

	done

done
wait
