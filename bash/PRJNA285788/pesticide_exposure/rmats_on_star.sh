#!/bin/bash

#SBATCH -A m2_jgu-funcpoll
#SBATCH -J rmats_star

#SBATCH -p devel
#SBATCH -n 3
#SBATCH -c 10
#SBATCH --mem-per-cpu 1G
#SBATCH -t 04:00:00

#SBATCH -v
#SBATCH -o rmats_star.%j.out
#SBATCH -e rmats_star.%j.err

module purge
module load bio/rmats-turbo/4.1.1-GCCcore-10.2.0

PROJECT_DIR=/lustre/project/m2_jgu-funcpoll
BAM_DIR=${PROJECT_DIR}/results/PRJNA285788/pesticide_exposure/after_trimming/alignment
RMATS_OUT_DIR=${PROJECT_DIR}/results/PRJNA285788/pesticide_exposure/after_trimming/diff_splicing
RUN_TABLE=${PROJECT_DIR}/data/datasets/PRJNA285788/accession/SraRunTable.csv
ANNOTATION_FILE=${PROJECT_DIR}/data/reference/annotation/GCA_907164935.1_iOsmBic2.1_genomic_liftoff.gtf

mkdir -p ${RMATS_OUT_DIR}/tmp_paths

EXP_GROUPS=("control" "thiacloprid" "imidacloprid")

LENGTH=${#EXP_GROUPS[@]}

for (( i=0; i<${LENGTH}-1; i++ ));
do
	for (( j=${i}+1; j<${LENGTH}; j++ ));
	do
		OUT_DIR=${RMATS_OUT_DIR}/${EXP_GROUPS[i]}_${EXP_GROUPS[j]}
		mkdir -p ${OUT_DIR}
		mkdir -p ${OUT_DIR}/tmp

		srun -n 1 -c 10 \
		python $EBROOTRMATSMINTURBO/rmats.py \
			--b1 ${RMATS_OUT_DIR}/tmp_paths/${EXP_GROUPS[i]}.txt \
			--b2 ${RMATS_OUT_DIR}/tmp_paths/${EXP_GROUPS[j]}.txt \
			--od ${OUT_DIR} --tmp ${OUT_DIR}/tmp \
			--gtf ${ANNOTATION_FILE} --readLength 100 \
			--variable-read-length --nthread 10 -t paired &
	done
done

wait
