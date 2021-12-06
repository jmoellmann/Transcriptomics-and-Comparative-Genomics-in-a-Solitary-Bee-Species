#!/bin/bash

#SBATCH -A m2_jgu-funcpoll
#SBATCH -J htseq_count

#SBATCH -p devel
#SBATCH -n 45
#SBATCH -c 3
#SBATCH -t 04:00:00
#SBATCH --mem-per-cpu 1G

#SBATCH -v
#SBATCH -o htseq_count.%j.out
#SBATCH -e htseq_count.%j.err

PROJECT_DIR=/lustre/project/m2_jgu-funcpoll
EXPERIMENT_DIR=${PROJECT_DIR}/data/datasets/PRJNA533756
DATA_DIR=${PROJECT_DIR}/results/PRJNA533756/after_trimming/alignment
ANNOTATION_FILE=${PROJECT_DIR}/data/reference/annotation/Apis_mellifera.Amel_HAv3.1.49.gtf

mkdir -p ${OUT_DIR}

module purge
module load bio/HTSeq/0.13.5-foss-2020a-Python-3.8.2

for GROUP in "adult_male" "virgin_queen" "mated_queen" "worker_forager" "worker_nurse";
do
	GROUP_ACCESSION=${EXPERIMENT_DIR}/accession/${GROUP}.txt
	GROUP_DATA_DIR=${DATA_DIR}/${GROUP}

	for SAMPLE in $(cat ${GROUP_ACCESSION});
	do
		if [ -f  ${GROUP_DATA_DIR}/${SAMPLE}_Log.final.out ]
		then 
			OUT_FILE=${GROUP_DATA_DIR}/${SAMPLE}_transcript_counts.txt

			touch ${OUT_FILE}

			BAM_FILE=${GROUP_DATA_DIR}/${SAMPLE}_Aligned.sortedByCoord.out.bam

			srun -n 1 -c 3 \
			htseq-count --nprocesses 3 --format bam --stranded no --order pos \
			--idattr transcript_id ${BAM_FILE} ${ANNOTATION_FILE} > ${OUT_FILE} &
		fi
	done

done
wait 
