#!/bin/bash

#SBATCH -A m2_jgu-funcpoll
#SBATCH -J fastqc_before_trimming

#SBATCH -p devel
#SBATCH -n 27
#SBATCH -c 1
#SBATCH -t 02:00:00
#SBATCH --mem-per-cpu 700M

#SBATCH -v
#SBATCH -o fastqc_before_trimming.%j.out
#SBATCH -e fastqc_before_trimming.%j.err

module purge
module load bio/FastQC/0.11.9-Java-11

PROJECT_DIR=/lustre/project/m2_jgu-funcpoll
EXPERIMENT_DIR=${PROJECT_DIR}/data/datasets/PRJEB9366

FASTQ_DATA_DIR=${EXPERIMENT_DIR}/fastq-files/
REPORT_OUT_DIR=${PROJECT_DIR}/results/PRJEB9366/before_trimming/quality/fastqc/

mkdir -p ${REPORT_OUT_DIR}

for FASTQ_FILE in $(ls ${FASTQ_DATA_DIR});
do
	SAMPLE=${FASTQ_FILE%".fastq.gz"}
	if [ ! -f ${REPORT_OUT_DIR}/${SAMPLE}_fastqc.html ];
	then

		IN_FILE=${FASTQ_DATA_DIR}/${FASTQ_FILE}

		srun -n 1 -c 1 --mem-per-cpu 700M \
		fastqc --outdir ${REPORT_OUT_DIR} ${IN_FILE} &
	fi
done
wait
