#!/bin/bash

#SBATCH -A m2_jgu-funcpoll
#SBATCH -J remove_bad_transcript


#SBATCH -p parallel
#SBATCH -n 24
#SBATCH -c 1
#SBATCH --mem-per-cpu 300M
#SBATCH -t 02:00:00

#SBATCH -v
#SBATCH -o remove_bad_transcript.%j.out
#SBATCH -e remove_bad_transcript.%j.err

module purge
module load bio/SAMtools/1.11-GCC-9.3.0    

PROJECT_DIR=/lustre/project/m2_jgu-funcpoll
EXPERIMENT_DIR=${PROJECT_DIR}/data/datasets/drone_vs_worker
TRANSCRIPT_TO_REMOVE=XR_003306190

for GROUP in drone worker;
do
	GROUP_ACCESSION=${EXPERIMENT_DIR}/accession/GSE85433_3pm_${GROUP}.txt
	BAM_FILE_DIR=${PROJECT_DIR}/results/drone_vs_worker/after_trimming/alignment/${GROUP}

	for SAMPLE in $(cat ${GROUP_ACCESSION});
	do
		
		BAM_IN_FILE=${BAM_FILE_DIR}/${SAMPLE}_Aligned.toTranscriptome.out.bam
		BAM_OUT_FILE=${BAM_FILE_DIR}/${SAMPLE}_Aligned.toTranscriptome.bad.transcripts.removed.out.bam
		
		srun -n 1 -c 1 \
		samtools view -h ${BAM_IN_FILE} | grep -v ${TRANSCRIPT_TO_REMOVE} | samtools view -bS -o ${BAM_OUT_FILE} &
		
	done

done
wait 
