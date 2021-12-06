#!/bin/bash

PROJECT_DIR=/lustre/project/m2_jgu-funcpoll
ENA_ACCESSION=${PROJECT_DIR}/data/datasets/PRJNA285788/accession/ena_accession.tsv

for GROUP in "male" "female";
do
	OUTPUT_PATH=${PROJECT_DIR}/data/datasets/PRJNA285788/fastq-files/${GROUP}
	ACCESSION_FILE=${PROJECT_DIR}/data/datasets/PRJNA285788/accession/${GROUP}.txt
	
	mkdir -p ${OUTPUT_PATH}

	for SAMPLE_ID in $(cat ${ACCESSION_FILE});
	do
		DOWNLOAD_LINK_1=$(cut -f 7 ${ENA_ACCESSION} | grep "${SAMPLE_ID}" | cut -f 1 -d ";")
		DOWNLOAD_LINK_2=$(cut -f 7 ${ENA_ACCESSION} | grep "${SAMPLE_ID}" | cut -f 2 -d ";")

		wget -N -c -o ${PROJECT_DIR}/log/${SAMPLE_ID}_1_log.txt -P ${OUTPUT_PATH} ${DOWNLOAD_LINK_1} &
		wget -N -c -o ${PROJECT_DIR}/log/${SAMPLE_ID}_2_log.txt -P ${OUTPUT_PATH} ${DOWNLOAD_LINK_2} &	
	done
done

wait

