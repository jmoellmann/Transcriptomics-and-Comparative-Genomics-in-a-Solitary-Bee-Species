#!/bin/bash

module purge

PROJECT_DIR=/lustre/project/m2_jgu-funcpoll/

for GROUP in "worker_nurse" "worker_forager" "adult_male" "virgin_queen" "mated_queen";
do
	ACCESSION_FILE=${PROJECT_DIR}/data/datasets/PRJNA533756/accession/${GROUP}.txt
	FASTQ_FILES=${PROJECT_DIR}/data/datasets/PRJNA533756/fastq-files/${GROUP}

	for SAMPLE_ID in $(cat ${ACCESSION_FILE});
	do
		if gzip -t ${FASTQ_FILES}/${SAMPLE_ID}.fastq.gz; then
			echo "${SAMPLE_ID} ok"
		else
			echo "${SAMPLE_ID} not ok"
		fi	
	done
done
