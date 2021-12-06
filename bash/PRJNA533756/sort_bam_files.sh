#!/bin/bash

#SBATCH -A m2_jgu-funcpoll
#SBATCH -J sort_bam_files


#SBATCH -p parallel
#SBATCH -n 48
#SBATCH -c 4
#SBATCH --mem-per-cpu 1G
#SBATCH -t 48:00:00

#SBATCH -v
#SBATCH -o sort_bam_files.%j.out
#SBATCH -e sort_bam_files.%j.err

module purge
module load bio/SAMtools/1.11-GCC-9.3.0 

PROJECT_DIR=/lustre/project/m2_jgu-funcpoll
EXPERIMENT_DIR=${PROJECT_DIR}/data/datasets/drone_vs_worker

for TIME_POINT in before_trimming after_trimming;
do

	for GROUP in drone worker;
	do
		GROUP_ACCESSION=${EXPERIMENT_DIR}/accession/GSE85433_3pm_${GROUP}.txt
		BAM_DATA_DIR=${PROJECT_DIR}/results/drone_vs_worker/${TIME_POINT}/alignment/${GROUP}

		for SAMPLE in $(cat ${GROUP_ACCESSION});
		do
			
			UNSORTED_BAM_FILE=${BAM_DATA_DIR}/${SAMPLE}_Aligned.out.bam
			SORTED_BAM_FILE=${BAM_DATA_DIR}/${SAMPLE}_Aligned_Sorted.out.bam.gz
			
			srun -n 1 -N 1 --exclusive -c 4 \
			samtools sort --threads 3 -m 1G ${UNSORTED_BAM_FILE} > ${SORTED_BAM_FILE} &
		done

	done
done
wait 
