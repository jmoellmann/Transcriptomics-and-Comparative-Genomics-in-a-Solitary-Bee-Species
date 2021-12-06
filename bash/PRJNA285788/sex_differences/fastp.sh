#!/bin/bash

#SBATCH -A m2_jgu-funcpoll
#SBATCH -J fastp

#SBATCH -p devel
#SBATCH -n 7
#SBATCH -c 1
#SBATCH -t 03:00:00
#SBATCH --mem-per-cpu 2G

#SBATCH -v
#SBATCH -o fastp.%j.out
#SBATCH -e fastp.%j.err

module purge
module load bio/fastp/0.20.1-GCC-8.3.0

PROJECT_DIR=/lustre/project/m2_jgu-funcpoll
EXPERIMENT_DIR=${PROJECT_DIR}/data/datasets/PRJNA285788
ADAPTERS=${PROJECT_DIR}/data/reference/adapters/illumina_adapters.fasta

for GROUP in "female" "male";
do
	GROUP_ACCESSION=${EXPERIMENT_DIR}/accession/${GROUP}.txt
	FASTQ_DATA_DIR=${EXPERIMENT_DIR}/fastq-files/${GROUP}
	FASTP_OUT_DIR=${EXPERIMENT_DIR}/fastp-trimmed/${GROUP}
	REPORT_OUT_DIR=${PROJECT_DIR}/results/PRJNA285788/after_trimming/quality/fastp/${GROUP}

	mkdir -p ${FASTP_OUT_DIR}
	mkdir -p ${REPORT_OUT_DIR}

	cd ${REPORT_OUT_DIR}

	for SAMPLE in $(cat ${GROUP_ACCESSION});
	do
		if [ ! -f ${REPORT_OUT_DIR}/${SAMPLE}_fastp.html ] || \
		[ ! -f  ${REPORT_OUT_DIR}/${SAMPLE}_fastp.json ];
		then

		srun -n 1 -c 1 --mem-per-cpu 2G \
		fastp --in1 ${FASTQ_DATA_DIR}/${SAMPLE}_1.fastq.gz \
		--in2 ${FASTQ_DATA_DIR}/${SAMPLE}_2.fastq.gz \
		--out1 ${FASTP_OUT_DIR}/${SAMPLE}_1_trimmed.fastq.gz \
		--out2 ${FASTP_OUT_DIR}/${SAMPLE}_2_trimmed.fastq.gz \
		--adapter_fasta ${ADAPTERS} \
		--length_required 50 \
		--json ${SAMPLE}_fastp.json \
		--html ${SAMPLE}_fastp.html \
		-w 1 &
		
		fi
	done

done
wait 
