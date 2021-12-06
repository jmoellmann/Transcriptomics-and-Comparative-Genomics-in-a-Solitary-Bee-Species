#!/bin/bash

#SBATCH -A m2_jgu-funcpoll
#SBATCH -J fastqc_before_trimming

#SBATCH -p devel
#SBATCH -n 7
#SBATCH -c 1
#SBATCH -t 02:00:00
#SBATCH --mem-per-cpu 700M

#SBATCH -v
#SBATCH -o fastqc_before_trimming.%j.out
#SBATCH -e fastqc_before_trimming.%j.err

module purge
module load bio/FastQC/0.11.9-Java-11 

PROJECT_DIR=/lustre/project/m2_jgu-funcpoll
EXPERIMENT_DIR=${PROJECT_DIR}/data/datasets/PRJNA285788

for GROUP in "female" "male";
do
	GROUP_ACCESSION=${EXPERIMENT_DIR}/accession/${GROUP}.txt
	DATA_DIR=${EXPERIMENT_DIR}/fastq-files/${GROUP}
	REPORT_OUT_DIR=${PROJECT_DIR}/results/PRJNA285788/before_trimming/quality/fastqc/${GROUP}

	mkdir -p ${REPORT_OUT_DIR}


	for SAMPLE in $(cat ${GROUP_ACCESSION});
	do
		if [ ! -f ${REPORT_OUT_DIR}/${SAMPLE}_fastqc.html ];
		then

			IN_FILE_1=${DATA_DIR}/${SAMPLE}_1.fastq.gz
			IN_FILE_2=${DATA_DIR}/${SAMPLE}_2.fastq.gz

			srun -n 1 -c 1 --mem-per-cpu 700M \
			fastqc --outdir ${REPORT_OUT_DIR} ${IN_FILE_1} &
			
			srun -n 1 -c 1 --mem-per-cpu 700M \
			fastqc --outdir ${REPORT_OUT_DIR} ${IN_FILE_2} &
		fi
	done

done
wait 
