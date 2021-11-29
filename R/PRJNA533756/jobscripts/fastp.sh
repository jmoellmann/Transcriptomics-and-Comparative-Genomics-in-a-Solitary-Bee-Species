#!/bin/bash

#SBATCH -A m2_jgu-funcpoll
#SBATCH -J fastp


#SBATCH -p smp
#SBATCH -n 18
#SBATCH -c 1
#SBATCH -C skylake
#SBATCH -t 03:00:00
#SBATCH --mem-per-cpu 2G

#SBATCH -v
#SBATCH -o fastp.%j.out
#SBATCH -e fastp.%j.err

module purge
module load bio/fastp/0.20.1-GCC-8.3.0

PROJECT_DIR=/lustre/project/m2_jgu-funcpoll
EXPERIMENT_DIR=${PROJECT_DIR}/data/datasets/PRJNA533756
ADAPTERS=${PROJECT_DIR}/data/reference/adapters/illumina_adapters.fasta

for GROUP in "worker_nurse" "worker_forager" "adult_male" "virgin_queen" "mated_queen";
do
	GROUP_ACCESSION=${EXPERIMENT_DIR}/accession/${GROUP}.txt
	FASTQ_DATA_DIR=${EXPERIMENT_DIR}/fastq-files/${GROUP}
	FASTP_OUT_DIR=${EXPERIMENT_DIR}/fastp-trimmed/${GROUP}
	REPORT_OUT_DIR=${PROJECT_DIR}/results/PRJNA533756/after_trimming/quality/fastp/${GROUP}

	mkdir -p ${FASTP_OUT_DIR}
	mkdir -p ${REPORT_OUT_DIR}

	cd ${REPORT_OUT_DIR}

	for SAMPLE in $(cat ${GROUP_ACCESSION});
	do
		if [ ! -f ${REPORT_OUT_DIR}/${SAMPLE}_fastp.html ] || \
		[ ! -f  ${REPORT_OUT_DIR}/${SAMPLE}_fastp.json ];
		then

		srun -n 1 -c 1 --mem-per-cpu 2G \
		fastp -i ${FASTQ_DATA_DIR}/${SAMPLE}.fastq.gz \
		-o ${FASTP_OUT_DIR}/${SAMPLE}_trimmed.fastq.gz \
		--adapter_fasta ${ADAPTERS} \
		--length_required 40 \
		--json ${SAMPLE}_fastp.json \
		--html ${SAMPLE}_fastp.html \
		-w 1 &
		
		fi
	done

done
wait 
