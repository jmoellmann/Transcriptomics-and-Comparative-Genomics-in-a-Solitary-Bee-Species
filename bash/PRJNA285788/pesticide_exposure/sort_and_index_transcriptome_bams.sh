#!/bin/bash

#SBATCH -A m2_jgu-funcpoll
#SBATCH -J htseq_count

#SBATCH -p devel
#SBATCH -n 45
#SBATCH -c 1
#SBATCH -t 04:00:00
#SBATCH --mem-per-cpu 1G

#SBATCH -v
#SBATCH -o htseq_count.%j.out
#SBATCH -e htseq_count.%j.err

module purge
module load bio/SAMtools/1.12-GCC-10.2.0

PROJECT_DIR=/lustre/project/m2_jgu-funcpoll
EXPERIMENT_DIR=${PROJECT_DIR}/data/datasets/PRJNA533756
DATA_DIR=${PROJECT_DIR}/results/PRJNA533756/after_trimming/alignment
ANNOTATION_FILE=${PROJECT_DIR}/data/reference/annotation/Apis_mellifera.Amel_HAv3.1.49.gtf

mkdir -p ${OUT_DIR}

for GROUP in "adult_male" "virgin_queen" "mated_queen" "worker_forager" "worker_nurse";
do
	GROUP_ACCESSION=${EXPERIMENT_DIR}/accession/${GROUP}.txt
	GROUP_DATA_DIR=${DATA_DIR}/${GROUP}

	for SAMPLE in $(cat ${GROUP_ACCESSION});
	do
		if [ -f  ${GROUP_DATA_DIR}/${SAMPLE}_Log.final.out ]
		then
			BAM_FILE=${GROUP_DATA_DIR}/${SAMPLE}_Aligned.toTranscriptome.out.bam
			SORTED_BAM_FILE=${GROUP_DATA_DIR}/${SAMPLE}_Aligned.toTranscriptome.sorted.out.bam

			srun -n 1 -c 1 samtools sort ${BAM_FILE} -o ${SORTED_BAM_FILE} \
			&& samtools index ${SORTED_BAM_FILE} &
		fi
	done

done
wait
