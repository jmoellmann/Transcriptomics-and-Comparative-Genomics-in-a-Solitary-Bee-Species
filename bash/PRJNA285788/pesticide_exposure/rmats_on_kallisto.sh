#!/bin/bash

#SBATCH -A m2_jgu-funcpoll
#SBATCH -J rmats_kallisto

#SBATCH -p devel
#SBATCH -n 45
#SBATCH -c 2
#SBATCH --mem-per-cpu 1G
#SBATCH -t 00:30:00

#SBATCH -v
#SBATCH -o rmats_kallisto.%j.out
#SBATCH -e rmats_kallisto.%j.err

module purge
module load bio/rmats-turbo/4.1.1-GCCcore-10.2.0

PROJECT_DIR=/lustre/project/m2_jgu-funcpoll
BAM_DIR=${PROJECT_DIR}/results/PRJNA533756/after_trimming/quantification/kallisto
RMATS_OUT_DIR=${PROJECT_DIR}/results/PRJNA533756/after_trimming/diff_splicing/kallisto
RUN_TABLE=${PROJECT_DIR}/data/datasets/PRJNA533756/accession/SraRunTable.csv
ANNOTATION_FILE=${PROJECT_DIR}/data/reference/annotation/Apis_mellifera.Amel_HAv3.1.49.gtf
BAM_NAME="pseudoalignments_with_nh_tag.bam"

mkdir -p ${RMATS_OUT_DIR}/tmp_paths

EXP_GROUPS=("adult_male" "virgin_queen" "mated_queen" "worker_forager" "worker_nurse")

LENGTH=${#EXP_GROUPS[@]}

for BODY_PART in "Thorax" "Abdomen" "Head"
do

	for EXP_GROUP in ${EXP_GROUPS[@]}
	do
		FILTER=($(cut -f 1,36,38 -d "," ${RUN_TABLE} | grep "${EXP_GROUP//_/ },${BODY_PART}" | cut -f 1 -d ","))
		SAMPLES=${RMATS_OUT_DIR}/tmp_paths/${BODY_PART}_${EXP_GROUP}.txt
		touch ${SAMPLES}
		echo -n "${BAM_DIR}/${FILTER[0]}/${BAM_NAME}","${BAM_DIR}/${FILTER[1]}/${BAM_NAME}","${BAM_DIR}/${FILTER[2]}/${BAM_NAME}" > ${SAMPLES}
	done

	for (( i=0; i<${LENGTH}-1; i++ ));
	do
		for (( j=${i}+1; j<${LENGTH}; j++ ));
		do
			OUT_DIR=${RMATS_OUT_DIR}/${BODY_PART}/${EXP_GROUPS[i]}_${EXP_GROUPS[j]}
			mkdir -p ${OUT_DIR}

			srun -n 1 -c 2 \
			python $EBROOTRMATSMINTURBO/rmats.py \
				--b1 ${RMATS_OUT_DIR}/tmp_paths/${BODY_PART}_${EXP_GROUPS[i]}.txt \
				--b2 ${RMATS_OUT_DIR}/tmp_paths/${BODY_PART}_${EXP_GROUPS[j]}.txt \
				--od ${OUT_DIR} --tmp ${OUT_DIR} \
				--gtf ${ANNOTATION_FILE} --readLength 50 \
				--variable-read-length --nthread 2 -t single &
		done
	done

done

wait
