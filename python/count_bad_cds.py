import os
from Bio import SeqIO

os.chdir("/media/jannik/Samsung_T5/masters_project/01_main_analysis/")

def count_bad_cds(cds_fasta):
   cds_seqs = [seq_rec for seq_rec in SeqIO.parse(cds_fasta, "fasta")]
   count = 0
   for cds_seq in cds_seqs:
      if len(cds_seq.seq) % 3 != 0:
         count += 1
   return count
   


cds_fasta_dir = "01_data/02_phylogenetic_analysis/fasta/01_cds_from_refseq/"

cds_fasta_file_paths = [cds_fasta_dir + species for species in os.listdir(cds_fasta_dir)]

for cds_fasta in cds_fasta_file_paths:
   print(count_bad_cds(cds_fasta), "bad cds fastas for species", os.path.split(cds_fasta)[-1][:-4])
