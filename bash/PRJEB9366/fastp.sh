#!/bin/bash

#SBATCH -A m2_jgu-funcpoll
#SBATCH -J fastp

#SBATCH -p devel
#SBATCH -n 27
#SBATCH -c 1
#SBATCH -t 03:00:00
#SBATCH --mem-per-cpu 2G

#SBATCH -v
#SBATCH -o fastp.%j.out
#SBATCH -e fastp.%j.err

module purge
module load bio/fastp/0.20.1-GCC-8.3.0

PROJECT_DIR=/lustre/project/m2_jgu-funcpoll
EXPERIMENT_DIR=${PROJECT_DIR}/data/datasets/PRJEB9366
ADAPTERS=${PROJECT_DIR}/data/reference/adapters/illumina_adapters.fasta

FASTQ_DATA_DIR=${EXPERIMENT_DIR}/fastq-files
FASTP_OUT_DIR=${EXPERIMENT_DIR}/fastp-trimmed
REPORT_OUT_DIR=${PROJECT_DIR}/results/PRJEB9366/after_trimming/quality/fastp

mkdir -p ${FASTP_OUT_DIR}
mkdir -p ${REPORT_OUT_DIR}

cd ${REPORT_OUT_DIR}

for FASTQ_FILE in $(ls ${FASTQ_DATA_DIR});
do
	SAMPLE=${FASTQ_FILE%".fastq.gz"}

	if [ ! -f ${REPORT_OUT_DIR}/${SAMPLE}_fastp.html ] || \
	[ ! -f  ${REPORT_OUT_DIR}/${SAMPLE}_fastp.json ];
	then

		srun -n 1 -c 1 \
		fastp -i ${FASTQ_DATA_DIR}/${FASTQ_FILE} \
		-o ${FASTP_OUT_DIR}/${SAMPLE}_trimmed.fastq.gz \
		--adapter_fasta ${ADAPTERS} \
		--length_required 40 \
		--json ${SAMPLE}_fastp.json \
		--html ${SAMPLE}_fastp.html \
		-w 1 &

	fi

done
wait
