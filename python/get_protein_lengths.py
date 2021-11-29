import os
import pysam
import pandas as pd

os.chdir("/media/jannik/Samsung_T5/masters_project/ApisTranscriptomics/")

def get_protein_length(fasta):
    fasta_file = pysam.FastaFile(fasta)
    df = pd.DataFrame()
    df["ProteinID"] = fasta_file.references
    df["Length"] = fasta_file.lengths
    return(df)

dm_fasta = "results/PRJNA285788/immune_gene_identification/primary_transcripts/Drosophila_melanogaster.faa"
ob_fasta = "results/PRJNA285788/immune_gene_identification/primary_transcripts/Osmia_bicornis.faa"
ol_fasta = "results/PRJNA285788/immune_gene_identification/primary_transcripts/Osmia_lignaria.faa"
bt_fasta = "results/PRJNA285788/immune_gene_identification/primary_transcripts/Bombus_terrestris.faa"
am_fasta = "results/PRJNA285788/immune_gene_identification/primary_transcripts/Apis_mellifera.faa"

get_protein_length(dm_fasta).to_csv("data/reference/annotation/dm_primary_transcript_lengths.csv")
get_protein_length(ob_fasta).to_csv("data/reference/annotation/ob_primary_transcript_lengths.csv")
get_protein_length(ol_fasta).to_csv("data/reference/annotation/ol_primary_transcript_lengths.csv")
get_protein_length(bt_fasta).to_csv("data/reference/annotation/bt_primary_transcript_lengths.csv")
get_protein_length(am_fasta).to_csv("data/reference/annotation/am_primary_transcript_lengths.csv")
