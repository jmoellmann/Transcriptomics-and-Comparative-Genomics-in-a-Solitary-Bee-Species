#!/bin/bash

#SBATCH -A m2_jgu-funcpoll
#SBATCH -J compress

#SBATCH -p smp
#SBATCH -n 2
#SBATCH -c 1
#SBATCH --mem-per-cpu 200M
#SBATCH -t 48:00:00

#SBATCH -v
#SBATCH -o compress.%j.out
#SBATCH -e compress.%j.err

module purge

PROJECT_DIR=/lustre/project/m2_jgu-funcpoll
DATA_DIR=${PROJECT_DIR}/data/datasets/GSE45408_pupa/fasterq-dump

for GROUP in drone worker queen;
do

	GROUP_DIR=${DATA_DIR}/${GROUP}

	for SAMPLE in $(ls ${GROUP_DIR});
	do
		if [[ ${SAMPLE} != *.gz ]];
		then
			srun -n 1 -c 1 \
			gzip ${GROUP_DIR}/${SAMPLE} &
		fi
	done

done
wait
