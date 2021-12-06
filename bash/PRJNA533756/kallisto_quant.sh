#!/bin/bash

#SBATCH -A m2_jgu-funcpoll
#SBATCH -J kallisto

#SBATCH -p devel
#SBATCH -n 3
#SBATCH -c 10
#SBATCH --mem-per-cpu 2G
#SBATCH -t 02:00:00

#SBATCH -v
#SBATCH -o kallisto.%j.out
#SBATCH -e kallisto.%j.err

module purge
module load bio/kallisto/0.46.1-foss-2019b

PROJECT_DIR=/lustre/project/m2_jgu-funcpoll
EXPERIMENT_DIR=${PROJECT_DIR}/data/datasets/PRJNA533756
INDEX_FILE=${PROJECT_DIR}/data/reference/kallisto_index/Apis_mellifera.49.kallisto_index.idx
KALLISTO_OUT_DIR=${PROJECT_DIR}/results/PRJNA533756/after_trimming/quantification/kallisto
ANNOTATION_FILE=${PROJECT_DIR}/data/reference/annotation/Apis_mellifera.Amel_HAv3.1.49.gtf
CHROMOSOMES_FILE=${PROJECT_DIR}/data/reference/annotation/chromosomes_and_lengths.txt

for GROUP in "adult_male" "worker_forager" "worker_nurse" "mated_queen" "virgin_queen";
do
	GROUP_ACCESSION=${EXPERIMENT_DIR}/accession/${GROUP}.txt
	FASTP_DATA_DIR=${EXPERIMENT_DIR}/fastp-trimmed/${GROUP}

	for SAMPLE in $(cat ${GROUP_ACCESSION});
	do
		if [ ! -f ${KALLISTO_OUT_DIR}/${SAMPLE}/abundance.h5 ];
		then
			mkdir -p ${KALLISTO_OUT_DIR}/${SAMPLE}

			IN_FILE=${FASTP_DATA_DIR}/${SAMPLE}_trimmed.fastq.gz

			srun -n 1 -c 10 \
			kallisto quant -i ${INDEX_FILE} -o ${KALLISTO_OUT_DIR}/${SAMPLE} --threads 10 --single -l 50 -s 1 ${IN_FILE} \
			--genomebam --gtf ${ANNOTATION_FILE} --chromosomes ${CHROMOSOMES_FILE} &> ${KALLISTO_OUT_DIR}/${SAMPLE}/${SAMPLE}_kallisto.log &
		fi

	done

done
wait
