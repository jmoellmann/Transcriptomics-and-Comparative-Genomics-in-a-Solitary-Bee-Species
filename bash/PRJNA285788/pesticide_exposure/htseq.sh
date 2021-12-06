#!/bin/bash

#SBATCH -A m2_jgu-funcpoll
#SBATCH -J htseq_count


#SBATCH -p smp
#SBATCH -n 24
#SBATCH -c 1
#SBATCH -t 04:00:00
#SBATCH --mem-per-cpu 2G

#SBATCH -v
#SBATCH -o htseq_count.%j.out
#SBATCH -e htseq_count.%j.err

module purge
module load bio/HTSeq/0.13.5-foss-2020a-Python-3.8.2

PROJECT_DIR=/lustre/project/m2_jgu-funcpoll
EXPERIMENT_DIR=${PROJECT_DIR}/data/datasets/drone_vs_worker
ANNOTATION_FILE=${PROJECT_DIR}/data/reference/annotation/Apis_mellifera.Amel_HAv3.1.49.gtf
OUT_DIR=${PROJECT_DIR}/results/drone_vs_worker/after_trimming/quantification/HTSeq

mkdir -p ${OUT_DIR}

for GROUP in drone worker;
do
	GROUP_ACCESSION=${EXPERIMENT_DIR}/accession/GSE85433_3pm_${GROUP}.txt
	BAM_DATA_DIR=${PROJECT_DIR}/results/drone_vs_worker/after_trimming/alignment/${GROUP}

	for SAMPLE in $(cat ${GROUP_ACCESSION});
	do
		if [ -f  ${BAM_DATA_DIR}/${SAMPLE}_Log.final.out ]
		then 
			
			OUT_FILE=${OUT_DIR}/${SAMPLE}_gene_counts.txt

			touch ${OUT_FILE}

			BAM_FILE=${BAM_DATA_DIR}/${SAMPLE}_Aligned.sortedByCoord.out.bam

			srun -n 1 -c 1 --mem-per-cpu 2G \
			htseq-count --format bam --order pos --stranded no ${BAM_FILE} ${ANNOTATION_FILE} > ${OUT_FILE} &
		fi
	done

done
wait 
