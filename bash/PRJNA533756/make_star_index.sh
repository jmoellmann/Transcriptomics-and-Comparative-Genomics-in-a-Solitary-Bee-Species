#!/bin/bash

#SBATCH -A m2_jgu-funcpoll
#SBATCH -J make_star_index


#SBATCH -p devel
#SBATCH -n 1
#SBATCH -c 10
#SBATCH --mem-per-cpu 2G

#SBATCH -v
#SBATCH -o make_star_index.%j.out
#SBATCH -e make_star_index.%j.err

module purge
module load bio/STAR/2.7.3a-GCC-9.3.0  

PROJECT_DIR=/lustre/project/m2_jgu-funcpoll
GENOME_DIR=${PROJECT_DIR}/data/reference/genome/ensembl/
OUT_DIR=${PROJECT_DIR}/data/reference/STAR_index
ANNOTATION_DIR=${PROJECT_DIR}/data/reference/annotation

mkdir -p ${OUT_DIR}

srun -n 1 -c 10 --exclusive \
STAR --runThreadN 10 --runMode genomeGenerate --genomeDir ${OUT_DIR} \
--genomeFastaFiles ${GENOME_DIR}/Apis_mellifera.Amel_HAv3.1.dna.primary_assembly.CM0099[31-47]*2.fa \
--sjdbGTFfile ${ANNOTATION_DIR}/Apis_mellifera.Amel_HAv3.1.49.gtf \
--sjdbOverhang 99 --genomeSAindexNbases 12

