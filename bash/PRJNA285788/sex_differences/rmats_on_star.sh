#!/bin/bash

#SBATCH -A m2_jgu-funcpoll
#SBATCH -J rmats_star

#SBATCH -p devel
#SBATCH -n 1
#SBATCH -c 10
#SBATCH --mem-per-cpu 1G
#SBATCH -t 04:00:00

#SBATCH -v
#SBATCH -o rmats_star.%j.out
#SBATCH -e rmats_star.%j.err

module purge
module load bio/rmats-turbo/4.1.1-GCCcore-10.2.0

PROJECT_DIR=/lustre/project/m2_jgu-funcpoll
BAM_DIR=${PROJECT_DIR}/results/PRJNA285788/sex_differences/after_trimming/alignment
RMATS_OUT_DIR=${PROJECT_DIR}/results/PRJNA285788/sex_differences/after_trimming/diff_splicing
RUN_TABLE=${PROJECT_DIR}/data/datasets/PRJNA285788/accession/SraRunTable.csv
ANNOTATION_FILE=${PROJECT_DIR}/data/reference/annotation/GCA_907164935.1_iOsmBic2.1_genomic_liftoff.gtf

mkdir -p ${RMATS_OUT_DIR}/tmp_paths

EXP_GROUPS=("male" "female")

srun -n 1 -c 10 \
python $EBROOTRMATSMINTURBO/rmats.py \
	--b1 ${RMATS_OUT_DIR}/tmp_paths/${EXP_GROUPS[0]}.txt \
	--b2 ${RMATS_OUT_DIR}/tmp_paths/${EXP_GROUPS[1]}.txt \
	--od ${RMATS_OUT_DIR} --tmp ${RMATS_OUT_DIR}/tmp_paths \
	--gtf ${ANNOTATION_FILE} --readLength 100 \
	--variable-read-length --nthread 10 -t paired &


wait
