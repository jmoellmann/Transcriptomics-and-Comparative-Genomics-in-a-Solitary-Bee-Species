#!/bin/bash

#SBATCH -A m2_jgu-funcpoll
#SBATCH -J rmats_star

#SBATCH -p parallel
#SBATCH -n 15
#SBATCH -c 1
#SBATCH --mem-per-cpu 1G
#SBATCH -t 00:30:00

#SBATCH -v
#SBATCH -o rmats_star.%j.out
#SBATCH -e rmats_star.%j.err

module purge
module load bio/rmats-turbo/4.1.1-GCCcore-10.2.0

PROJECT_DIR=/lustre/project/m2_jgu-funcpoll
BAM_DIR=${PROJECT_DIR}/results/PRJNA533756/after_trimming/alignment
RMATS_OUT_DIR=${PROJECT_DIR}/results/PRJNA533756/after_trimming/diff_splicing/STAR
RUN_TABLE=${PROJECT_DIR}/data/datasets/PRJNA533756/accession/SraRunTable.csv
ANNOTATION_FILE=${PROJECT_DIR}/data/reference/annotation/Apis_mellifera.Amel_HAv3.1.49.gtf

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
		find ${BAM_DIR}/${EXP_GROUP}/*.bam | grep -e "${FILTER[0]}" -e "${FILTER[1]}" -e "${FILTER[2]}" > ${SAMPLES}
		sed -z -i 's/\n//3; s/\n/,/g' ${SAMPLES}
	done

	for (( i=0; i<${LENGTH}-1; i++ ));
	do
		for (( j=${i}+1; j<${LENGTH}; j++ ));
		do
			OUT_DIR=${RMATS_OUT_DIR}/${BODY_PART}/${EXP_GROUPS[i]}_${EXP_GROUPS[j]}
			mkdir -p ${OUT_DIR}

			srun -n 1 -c 1 \
			python $EBROOTRMATSMINTURBO/rmats.py \
				--b1 ${RMATS_OUT_DIR}/tmp_paths/${BODY_PART}_${EXP_GROUPS[i]}.txt \
				--b2 ${RMATS_OUT_DIR}/tmp_paths/${BODY_PART}_${EXP_GROUPS[j]}.txt \
				--od ${OUT_DIR} --tmp ${OUT_DIR} \
				--gtf ${ANNOTATION_FILE} --readLength 50 \
				--variable-read-length --nthread 1 -t single &
		done
	done

done

wait
