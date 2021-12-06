#!/bin/bash

#SBATCH -A m2_jgu-funcpoll
#SBATCH -J qualimap_after_trimming


#SBATCH -p parallel
#SBATCH -n 24
#SBATCH -c 4
#SBATCH -t 08:00:00
#SBATCH --mem-per-cpu 2G

#SBATCH -v
#SBATCH -o qualimap_after_trimming.%j.out
#SBATCH -e qualimap_after_trimming.%j.err

module purge
module load bio/Qualimap/2.2.1-foss-2019a-R-3.6.0

PROJECT_DIR=/lustre/project/m2_jgu-funcpoll
EXPERIMENT_DIR=${PROJECT_DIR}/data/datasets/drone_vs_worker
ANNOTATION_FILE=${PROJECT_DIR}/data/reference/annotation/Apis_mellifera.Amel_HAv3.1.49.gtf

for GROUP in drone worker;
do
	GROUP_ACCESSION=${EXPERIMENT_DIR}/accession/GSE85433_3pm_${GROUP}.txt
	BAM_DATA_DIR=${PROJECT_DIR}/results/drone_vs_worker/after_trimming/alignment/${GROUP}
	BAM_QC_OUT_DIR=${PROJECT_DIR}/results/drone_vs_worker/after_trimming/quality/qualimap/bamqc/${GROUP}
	RNA_SEQ_QC_OUT_DIR=${PROJECT_DIR}/results/drone_vs_worker/after_trimming/quality/qualimap/rnaseq/${GROUP}

	for SAMPLE in $(cat ${GROUP_ACCESSION});
	do

		mkdir -p ${BAM_QC_OUT_DIR}/${SAMPLE}
		mkdir -p ${RNA_SEQ_QC_OUT_DIR}/${SAMPLE}

		BAM_FILE=${BAM_DATA_DIR}/${SAMPLE}_Aligned_Sorted.out.bam.gz

		srun -n 1 -N 1 -c 4 --mem-per-cpu 2G \
		qualimap bamqc -bam ${BAM_FILE} -gff ${ANNOTATION_FILE} \
		-outdir ${BAM_QC_OUT_DIR}/${SAMPLE} -c \
		--java-mem-size=7G &

		: '
		srun -n 1 -N 1 --exclusive -c 2 --mem-per-cpu 2G \
		qualimap rnaseq -bam ${BAM_FILE} -gtf ${ANNOTATION_FILE} \
		-outdir ${RNA_SEQ_QC_OUT_DIR}/${SAMPLE} -oc ${SAMPLE}_gene_counts.txt -pe \
		--java-mem-size=4G &
		'		
		
	done

done
wait 
