#!/bin/bash

PROJECT_DIR=/media/jannik/Samsung_T5/masters_project/01_main_analysis
PHYLO_DIR=${PROJECT_DIR}/01_data/02_phylogenetic_analysis
CDS_ALIGNMENTS_DIR=${PHYLO_DIR}/fasta/06_cds_alignments_v2

mkdir -p ${CDS_ALIGNMENTS_DIR}
touch ./log.txt
touch ./progress.txt

for ORTHOGROUP in $(ls ${PHYLO_DIR}/fasta/04_sco_protein_alignments_v2/);
do
	echo ${ORTHOGROUP} >> progress.txt
	perl ${PROJECT_DIR}/03_src/pal2nal.v14/pal2nal.pl ${PHYLO_DIR}/fasta/04_sco_protein_alignments_v2/${ORTHOGROUP} ${PHYLO_DIR}/fasta/05_cds_for_alignments_v2/${ORTHOGROUP} -output fasta > ${CDS_ALIGNMENTS_DIR}/${ORTHOGROUP} 2>> ./log.txt &
done

wait
