#! /bin/bash

PROJECT_DIR=/media/jannik/Samsung_T5/masters_project/01_main_analysis
SMASHPP=${PROJECT_DIR}/03_src/smashpp/experiment/bin/smashpp
CHROMOSOMES_DIR=${PROJECT_DIR}/01_data/02_phylogenetic_analysis/chr_rearrangements/chrs
TMP_DIR=${PROJECT_DIR}/01_data/02_phylogenetic_analysis/chr_rearrangements/tmp
OUT_DIR=${PROJECT_DIR}/01_data/02_phylogenetic_analysis/chr_rearrangements/svg
CWD=$(pwd)

mkdir -p ${OUT_DIR}
mkdir -p ${TMP_DIR}

cd ${TMP_DIR}

CHROMOSOMES=(${CHROMOSOMES_DIR}/*)

LENGTH=${#CHROMOSOMES[@]}

for (( i = 0; i < ${LENGTH}-1; i++ ));
do
	for (( j = ${i}+1; j < ${LENGTH}; j++ ));
	do
		SEQ_1=${CHROMOSOMES[i]}
		SEQ_2=${CHROMOSOMES[j]}
		${SMASHPP} -r ${SEQ_1} -t ${SEQ_2} -m 1000 -f 10 -th 1.2 -l 0 -n 2 -sb yes -ar yes
		#${SMASHPP} -viz -o ${OUT_DIR}/${SEQ_1##*/}.${SEQ_2##*/}.svg ${TMP_DIR}/${SEQ_1##*/}.${SEQ_2##*/}.pos
	done
done

cd ${CWD}
