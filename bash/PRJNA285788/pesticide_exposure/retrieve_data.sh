#!/bin/bash

#SBATCH -A m2_jgu-funcpoll
#SBATCH -J retrieve_data

#SBATCH -p smp
#SBATCH -n 8
#SBATCH -c 1
#SBATCH --mem-per-cpu 2G
#SBATCH -t 48:00:00

#SBATCH -v
#SBATCH -o retrieve_data.%j.out
#SBATCH -e retrieve_data.%j.err

module purge
module load bio/SRA-Toolkit/2.10.2-centos_linux64

PROJECT_DIR=/lustre/project/m2_jgu-funcpoll/

for GROUP in worker drone queen;
do
	OUTPUT_PATH=${PROJECT_DIR}/data/datasets/GSE45408_pupa/fasterq-dump/${GROUP}
	ACCESSION_FILE=${PROJECT_DIR}/data/datasets/GSE45408_pupa/accession/${GROUP}.txt
	mkdir -p ${OUTPUT_PATH}

	for SAMPLE_ID in $(cat ${ACCESSION_FILE});
	do
		if [ ! -f ${OUTPUT_PATH}/${SAMPLE_ID}.fastq.gz ]
		then
			echo ${SAMPLE_ID} && \
			srun -n 1 -c 1 \
			fastq-dump --gzip -O ${OUTPUT_PATH} ${SAMPLE_ID} &
		fi
	done
done

module purge
