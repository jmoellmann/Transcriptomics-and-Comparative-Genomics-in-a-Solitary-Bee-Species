#!/bin/bash

module purge
module load bio/SRA-Toolkit/2.10.2-centos_linux64

PROJECT_DIR=/lustre/project/m2_jgu-funcpoll
ENA_ACCESSION=${PROJECT_DIR}/data/datasets/PRJNA533756/accession/ena_accession.tsv

for GROUP in "worker_nurse" "worker_forager" "adult_male" "virgin_queen" "mated_queen";
do
	OUTPUT_PATH=${PROJECT_DIR}/data/datasets/PRJNA533756/fastq-files/${GROUP}
	ACCESSION_FILE=${PROJECT_DIR}/data/datasets/PRJNA533756/accession/${GROUP}.txt
	
	mkdir -p ${OUTPUT_PATH}

	for SAMPLE_ID in $(cat ${ACCESSION_FILE});
	do
		if [ ! -f ${OUTPUT_PATH}/${SAMPLE_ID}.fastq.gz ]
		then
			DOWNLOAD_LINK=$(cut -f 7 ${ENA_ACCESSION} | grep "${SAMPLE_ID}")			
			wget -P ${OUTPUT_PATH} ${DOWNLOAD_LINK} &
		fi	
	done
done

wait

module purge
