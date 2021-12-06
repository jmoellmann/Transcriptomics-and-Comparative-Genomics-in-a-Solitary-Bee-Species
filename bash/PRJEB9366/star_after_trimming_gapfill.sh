#!/bin/bash

#SBATCH -A m2_jgu-funcpoll
#SBATCH -J star_after_trimming

#SBATCH -C skylake
#SBATCH -p smp
#SBATCH -n 1
#SBATCH -c 2
#SBATCH --mem-per-cpu 2G
#SBATCH -t 12:00:00

#SBATCH -v
#SBATCH -o star_after_trimming.%j.out
#SBATCH -e star_after_trimming.%j.err

module purge
module load bio/STAR/2.7.3a-GCC-9.3.0  

PROJECT_DIR=/lustre/project/m2_jgu-funcpoll
EXPERIMENT_DIR=${PROJECT_DIR}/data/datasets/GSE145395_queen
INDEX_DIR=${PROJECT_DIR}/data/reference/STAR_index

for GROUP in queen;
do
	GROUP_ACCESSION=${EXPERIMENT_DIR}/accession/${GROUP}.txt
	FASTP_DATA_DIR=${EXPERIMENT_DIR}/fastp-trimmed/${GROUP}
	STAR_OUT_DIR=${PROJECT_DIR}/results/GSE145395_queen/after_trimming/alignment/${GROUP}

	mkdir -p ${STAR_OUT_DIR}

	for SAMPLE in $(cat ${GROUP_ACCESSION});
	do

		if [ ! -f  ${STAR_OUT_DIR}/${SAMPLE}_Log.out ]
		then 

			#echo ${SAMPLE}		

			IN_FILE=${FASTP_DATA_DIR}/${SAMPLE}_trimmed.fastq
			
			srun -n 1 -c 2 \
			STAR --runThreadN 2 --genomeDir ${INDEX_DIR} \
			--readFilesIn ${IN_FILE} \
			--genomeLoad LoadAndRemove --outFileNamePrefix ${STAR_OUT_DIR}/${SAMPLE}_ \
			--outSAMtype BAM SortedByCoordinate --quantMode GeneCounts \
			--limitBAMsortRAM 17000000000 &

		fi
		
	done

done
wait 
