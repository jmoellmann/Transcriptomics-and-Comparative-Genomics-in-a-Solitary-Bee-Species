import re
import os
import HTSeq
import pandas as pd
from tqdm import tqdm

os.chdir("/media/jannik/Samsung_T5/masters_project/01_main_analysis/")

def rename_cds(cds_alignment_path, cds_path):
   out_dir = "01_data/02_phylogenetic_analysis/fasta/06_cds_alignments_renamed_v2/"
   orthogroup_accession = os.path.split(cds_alignment_path)[-1]
   out_fasta_file = open(out_dir + orthogroup_accession, "w")
   aligned_sequences = [seq for seq in HTSeq.FastaReader(cds_alignment_path)]
   unaligned_sequences = [seq for seq in HTSeq.FastaReader(cds_path)]
   for aligned_sequence in aligned_sequences:
      for unaligned_sequence in unaligned_sequences:
         if bool(re.search(aligned_sequence.name, unaligned_sequence.name)):
            break
      species_name = unaligned_sequence.name.split("|")[0][:-5]
      aligned_sequence.name = species_name + "_" + aligned_sequence.name
      aligned_sequence.write_to_fasta_file(out_fasta_file)
      #except Exception:
         #out_fasta_file.write("ERROR\n")
   out_fasta_file.close()
   
aligned_cds_dir = "01_data/02_phylogenetic_analysis/fasta/06_cds_alignments_v2/"
unaligned_cds_dir = "01_data/02_phylogenetic_analysis/fasta/05_cds_for_alignments_v2/"

aligned_cds_file_paths = [aligned_cds_dir + i for i in os.listdir(aligned_cds_dir)]
unaligned_cds_file_paths = [unaligned_cds_dir + i for i in os.listdir(unaligned_cds_dir)]

try:
   os.mkdir("01_data/02_phylogenetic_analysis/fasta/06_cds_alignments_renamed_v2/")
except:
   None
   
for i in tqdm(range(len(aligned_cds_file_paths)), desc = "Processing MSAs:"):
   cds_alignment_path = aligned_cds_file_paths[i]
   cds_path = unaligned_cds_file_paths[i]
   rename_cds(cds_alignment_path, cds_path)

