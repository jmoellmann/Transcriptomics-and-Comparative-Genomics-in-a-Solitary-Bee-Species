#! usr/bin/env bash

PROJECT_DIR=/lustre/project/m2_jgu-funcpoll
FILES_FOLDER=${PROJECT_DIR}/results/drone_vs_worker/after_trimming/quantification/HTSeq


for FILE in $(ls ${FILES_FOLDER});
do

	OUT_FILE=$(echo ${FILE} | cut -f 1 -d .)_prefix_removed.txt

	sed 's/gene://g' ${FILES_FOLDER}/${FILE} > ${FILES_FOLDER}/${OUT_FILE}
done

