#!/bin/bash

#SBATCH -A m2_jgu-funcpoll
#SBATCH -J star_after_trimming

#SBATCH -p devel
#SBATCH -n 7
#SBATCH -c 5
#SBATCH --mem-per-cpu 2G
#SBATCH -t 04:00:00

#SBATCH -v
#SBATCH -o star_after_trimming.%j.out
#SBATCH -e star_after_trimming.%j.err

module purge
module load bio/STAR/2.7.3a-GCC-9.3.0  

PROJECT_DIR=/lustre/project/m2_jgu-funcpoll
EXPERIMENT_DIR=${PROJECT_DIR}/data/datasets/PRJNA285788
INDEX_DIR=${PROJECT_DIR}/data/reference/STAR_index_iOsmBic2.1

for GROUP in "female" "male";
do
	GROUP_ACCESSION=${EXPERIMENT_DIR}/accession/${GROUP}.txt
	DATA_DIR=${EXPERIMENT_DIR}/sex_differences/fastp-trimmed/${GROUP}
	STAR_OUT_DIR=${PROJECT_DIR}/results/PRJNA285788/sex_differences/after_trimming/alignment/${GROUP}

	mkdir -p ${STAR_OUT_DIR}

	for SAMPLE in $(cat ${GROUP_ACCESSION});
	do

		if [ ! -f  ${STAR_OUT_DIR}/${SAMPLE}_Log.final.out ]
		then 
		
			IN_FILE_1=${DATA_DIR}/${SAMPLE}_1_trimmed.fastq.gz
			IN_FILE_2=${DATA_DIR}/${SAMPLE}_2_trimmed.fastq.gz			

			srun -n 1 -c 5 \
			STAR --runThreadN 5 --genomeDir ${INDEX_DIR} \
			--readFilesIn ${IN_FILE_1} ${IN_FILE_2} --readFilesCommand gunzip -c \
			--genomeLoad NoSharedMemory --outFileNamePrefix ${STAR_OUT_DIR}/${SAMPLE}_ \
			--outSAMtype BAM Unsorted \
			--limitBAMsortRAM 8500000000 &

		fi
		
	done

done
wait 
