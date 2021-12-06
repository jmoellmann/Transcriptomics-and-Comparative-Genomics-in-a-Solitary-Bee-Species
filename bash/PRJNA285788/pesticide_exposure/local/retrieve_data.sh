#!/bin/bash

PROJECT_DIR=/lustre/project/m2_jgu-funcpoll
ENA_ACCESSION=${PROJECT_DIR}/data/datasets/PRJNA285788/accession/ena_accession.tsv

for GROUP in "control" "imidacloprid" "thiacloprid";
do
	OUTPUT_PATH=${PROJECT_DIR}/data/datasets/PRJNA285788/pesticide_exposure/fastq-files/${GROUP}
	ACCESSION_FILE=${PROJECT_DIR}/data/datasets/PRJNA285788/accession/${GROUP}.txt
	
	mkdir -p ${OUTPUT_PATH}

	for SAMPLE_ID in $(cat ${ACCESSION_FILE});
	do
		DOWNLOAD_LINK_1=$(cut -f 7 ${ENA_ACCESSION} | grep "${SAMPLE_ID}" | cut -f 1 -d ";")
		DOWNLOAD_LINK_2=$(cut -f 7 ${ENA_ACCESSION} | grep "${SAMPLE_ID}" | cut -f 2 -d ";")

		curl ${DOWNLOAD_LINK_1} > ${OUTPUT_PATH}/${SAMPLE_ID}_1.fastq.gz 2> ${PROJECT_DIR}/log/${SAMPLE_ID}_1_log.txt &
		curl ${DOWNLOAD_LINK_2} > ${OUTPUT_PATH}/${SAMPLE_ID}_2fastq.gz 2> ${PROJECT_DIR}/log/${SAMPLE_ID}_2_log.txt &	
	done
done

wait

