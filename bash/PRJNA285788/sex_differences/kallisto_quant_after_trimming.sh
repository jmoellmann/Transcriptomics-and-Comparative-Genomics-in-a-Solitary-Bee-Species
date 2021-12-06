#!/bin/bash

#SBATCH -A m2_jgu-funcpoll
#SBATCH -J kallisto

#SBATCH -p devel
#SBATCH -n 7
#SBATCH -c 5
#SBATCH --mem-per-cpu 2G
#SBATCH -t 02:00:00

#SBATCH -v
#SBATCH -o kallisto.%j.out
#SBATCH -e kallisto.%j.err

module purge
module load bio/kallisto/0.46.1-foss-2019b

PROJECT_DIR=/lustre/project/m2_jgu-funcpoll
EXPERIMENT_DIR=${PROJECT_DIR}/data/datasets/PRJNA285788
INDEX_FILE=${PROJECT_DIR}/data/reference/kallisto_index/Osmia_bicornis_v3.kallisto_index.idx
OUT_DIR=${PROJECT_DIR}/results/PRJNA285788/after_trimming/quantification/kallisto
ANNOTATION_FILE=${PROJECT_DIR}/data/reference/annotation/GCF_004153925.1_Obicornis_v3_genomic.gtf

for GROUP in "female" "male";
do
	GROUP_ACCESSION=${EXPERIMENT_DIR}/accession/${GROUP}.txt
	DATA_DIR=${EXPERIMENT_DIR}/fastp-trimmed/${GROUP}

	for SAMPLE in $(cat ${GROUP_ACCESSION});
	do
		if [ ! -f ${OUT_DIR}/${SAMPLE}/abundance.h5 ];
		then
			mkdir -p ${OUT_DIR}/${SAMPLE}

			IN_FILE_1=${DATA_DIR}/${SAMPLE}_1_trimmed.fastq.gz
			IN_FILE_2=${DATA_DIR}/${SAMPLE}_2_trimmed.fastq.gz

			srun -n 1 -c 5 \
			kallisto quant -i ${INDEX_FILE} -o ${OUT_DIR}/${SAMPLE} --threads 5 ${IN_FILE_1} \
			${IN_FILE_2} &> ${OUT_DIR}/${SAMPLE}/${SAMPLE}_kallisto.log &
		fi

	done

done
wait
