#!/bin/bash

#SBATCH -A m2_jgu-funcpoll
#SBATCH -J annotate_bams

#SBATCH -p devel
#SBATCH -n 45
#SBATCH -c 1
#SBATCH --mem-per-cpu 2G
#SBATCH -t 00:30:00

#SBATCH -v
#SBATCH -o annotate_bams.%j.out
#SBATCH -e annotate_bams.%j.err

module purge
module load bio/SAMtools/1.12-GCC-10.2.0

PROJECT_DIR=/lustre/project/m2_jgu-funcpoll
BAM_DIR=${PROJECT_DIR}/results/PRJNA533756/after_trimming/quantification/kallisto

for SAMPLE in $(ls ${BAM_DIR});
do
	srun -n 1 -c 1 \
	samtools view -h ${BAM_DIR}/${SAMPLE}/pseudoalignments.bam | awk -v OFS='\t' '{print $0, "NH:i:1"}' | \
	samtools view -b > ${BAM_DIR}/${SAMPLE}/pseudoalignments_with_nh_tag.bam &
done

wait
module purge

