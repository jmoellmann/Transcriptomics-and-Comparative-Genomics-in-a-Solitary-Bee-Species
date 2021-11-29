#!/bin/bash

#SBATCH -A m2_jgu-funcpoll
#SBATCH -J fastqc_after_trimming


#SBATCH -p smp
#SBATCH -n 18
#SBATCH -c 1
#SBATCH -C skylake
#SBATCH -t 02:00:00
#SBATCH --mem-per-cpu 700M

#SBATCH -v
#SBATCH -o fastqc_after_trimming.%j.out
#SBATCH -e fastqc_after_trimming.%j.err

module purge
module load bio/FastQC/0.11.9-Java-11

PROJECT_DIR=/lustre/project/m2_jgu-funcpoll
EXPERIMENT_DIR=${PROJECT_DIR}/data/datasets/PRJNA533756

for GROUP in "worker_nurse" "worker_forager" "adult_male" "virgin_queen" "mated_queen";
do
	GROUP_ACCESSION=${EXPERIMENT_DIR}/accession/${GROUP}.txt
	FASTQ_DATA_DIR=${EXPERIMENT_DIR}/fastp-trimmed/${GROUP}
	REPORT_OUT_DIR=${PROJECT_DIR}/results/PRJNA533756/after_trimming/quality/fastqc/${GROUP}

	mkdir -p ${REPORT_OUT_DIR}


	for SAMPLE in $(cat ${GROUP_ACCESSION});
	do
		if [ ! -f ${REPORT_OUT_DIR}/${SAMPLE}_fastqc.html ];
		then

			IN_FILE=${FASTQ_DATA_DIR}/${SAMPLE}_trimmed.fastq.gz

			srun -n 1 -c 1 --mem-per-cpu 700M \
			fastqc --outdir ${REPORT_OUT_DIR} ${IN_FILE} &
		fi	
	done

done
wait
