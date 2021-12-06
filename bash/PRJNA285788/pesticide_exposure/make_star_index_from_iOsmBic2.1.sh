#!/bin/bash

#SBATCH -A m2_jgu-funcpoll
#SBATCH -J make_star_index

#SBATCH -p devel
#SBATCH -n 1
#SBATCH -c 20
#SBATCH --mem-per-cpu 2G

#SBATCH -v
#SBATCH -o make_star_index.%j.out
#SBATCH -e make_star_index.%j.err

module purge
module load bio/STAR/2.7.3a-GCC-9.3.0  

PROJECT_DIR=/lustre/project/m2_jgu-funcpoll
GENOME_DIR=${PROJECT_DIR}/data/reference/genome/ncbi/
ANNOTATION_DIR=${PROJECT_DIR}/data/reference/annotation
OUT_DIR=${PROJECT_DIR}/data/reference/STAR_index_iOsmBic2.1

mkdir -p ${OUT_DIR}

srun -n 1 -c 20 --exclusive \
STAR --runThreadN 20 --runMode genomeGenerate --genomeDir ${OUT_DIR} \
--genomeFastaFiles ${GENOME_DIR}/GCA_907164935.1_iOsmBic2.1_genomic.fna \
--sjdbGTFfile ${ANNOTATION_DIR}/GCA_907164935.1_iOsmBic2.1_genomic_liftoff.gtf \
--sjdbOverhang 99 --genomeSAindexNbases 12 --genomeChrBinNbits 14

