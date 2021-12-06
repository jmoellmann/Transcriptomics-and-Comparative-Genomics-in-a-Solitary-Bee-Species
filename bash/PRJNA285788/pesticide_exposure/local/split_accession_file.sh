#!/bin/bash

module purge

PROJECT_DIR=/lustre/project/m2_jgu-funcpoll

METADATA=${PROJECT_DIR}/data/datasets/PRJNA533756/accession/SraRunTable.csv

for GROUP in "worker nurse" "worker forager" "adult male" "virgin queen" "mated queen";
do

	GROUP_ID=${GROUP//" "/"_"}
	GROUP_ACCESSION_FILE=${PROJECT_DIR}/data/datasets/PRJNA533756/accession/${GROUP_ID}.txt
	
	if [ ! -f ${GROUP_ACCESSION_FILE} ]; then
		cut -f 1,36 -d , ${METADATA} | grep "${GROUP}" | cut -f 1 -d , > ${GROUP_ACCESSION_FILE}
	fi
done

