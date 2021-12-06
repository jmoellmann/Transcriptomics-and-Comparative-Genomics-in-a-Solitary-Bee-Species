#!/bin/bash

#SBATCH -A m2_jgu-funcpoll
#SBATCH -J compress

#SBATCH -p parallel
#SBATCH -n 9
#SBATCH -c 1
#SBATCH --mem-per-cpu 200M
#SBATCH -t 02:00:00

#SBATCH -v
#SBATCH -o compress.%j.out
#SBATCH -e compress.%j.err

module purge

module load

PROJECT_DIR=/lustre/project/m2_jgu-funcpoll
DATA_DIR=${PROJECT_DIR}/data/datasets/GSE45408_pupa/fasterq-dump

for GROUP in drone worker queen;
do

	GROUP_DIR=${DATA_DIR}/${GROUP}

	for SAMPLE in $(ls ${GROUP_DIR});
	do
		if [[ ${SAMPLE} != *.gz ]];
		then
			srun 
			gzip ${GROUP_DIR}/${SAMPLE} &
		fi
	done

done
wait
