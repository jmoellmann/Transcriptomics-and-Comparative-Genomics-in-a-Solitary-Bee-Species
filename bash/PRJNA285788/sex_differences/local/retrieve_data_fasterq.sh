#!/bin/bash

module purge
module load bio/SRA-Toolkit/2.10.2-centos_linux64

PROJECT_DIR=/lustre/project/m2_jgu-funcpoll

for GROUP in "male" "female";
do
	OUTPUT_PATH=${PROJECT_DIR}/data/datasets/PRJNA285788/fasterq-dump/${GROUP}
	ACCESSION_FILE=data/datasets/PRJNA285788/accession/${GROUP}.txt
	
	mkdir -p ${OUTPUT_PATH}

	for SAMPLE_ID in $(cat ${ACCESSION_FILE});
	do
		fasterq-dump -p -O ${OUTPUT_PATH} ${SAMPLE_ID} &
	done
done

wait
module purge
