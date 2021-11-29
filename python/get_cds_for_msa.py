import os
import re
import HTSeq
#qimport pandas as pd
from tqdm import tqdm

os.chdir("/media/jannik/Samsung_T5/masters_project/01_main_analysis/")

def get_cds_for_msa(msa_path, coding_sequences):
   out_dir = "01_data/02_phylogenetic_analysis/fasta/05_cds_for_alignments_v2/"
   out_fasta_file = open(out_dir + os.path.split(msa_path)[-1], "w")
   sequence_ids = [seq.name for seq in HTSeq.FastaReader(msa_path)]
   for sequence_id in sequence_ids:
      try:
         coding_sequences[sequence_id].write_to_fasta_file(out_fasta_file)
      except Exception:
         out_fasta_file.write("ERROR\n")
   out_fasta_file.close()
   
def get_cds_for_msa_with_mrna(msa_path, coding_sequences):
   out_dir = "01_data/02_phylogenetic_analysis/fasta/05_cds_for_alignments_v2/"
   out_fasta_file = open(out_dir + os.path.split(msa_path)[-1], "w")
   sequence_ids = [seq.name for seq in HTSeq.FastaReader(msa_path)]
   for sequence_id in sequence_ids:
      accession = prot2transcript[sequence_id][1]
      try:
         coding_sequences[accession].write_to_fasta_file(out_fasta_file)
      except Exception:
         out_fasta_file.write("ERROR\n")
   out_fasta_file.close()
   
def gather_cds_files(cds_file_paths):
   coding_sequences = dict()
   for cds_file_path in cds_file_paths:
      species_name = os.path.split(cds_file_path)[-1][:-4]
      temp_dict_v1 = dict( (seq.name, seq) for seq in HTSeq.FastaReader(cds_file_path) )
      temp_dict_v2 = dict()
      for key, value in temp_dict_v1.items():
         temp_key = re.search("[XNY]P_[0-9]+\\.[0-9]", key)[0]
         temp_seq = value
         temp_seq.name = species_name + "_" + temp_key
         temp_dict_v2[temp_key] = temp_seq
      coding_sequences.update(temp_dict_v2)
   return coding_sequences

def gather_cds_files_from_mrna(cds_file_paths):
   coding_sequences = dict()
   for cds_file_path in cds_file_paths:
      species_name = os.path.split(cds_file_path)[-1][:-4]
      temp_dict = dict( (seq.name, seq) for seq in HTSeq.FastaReader(cds_file_path) )
      for key, value in temp_dict.items():
         temp_seq = value
         temp_seq.name = species_name + "_" + key
         temp_dict[key] = temp_seq
      coding_sequences.update(temp_dict)
   return coding_sequences


cds_dir = "01_data/02_phylogenetic_analysis/fasta/01_cds_from_refseq/"
cds_file_paths = [cds_dir + i for i in os.listdir(cds_dir)]
coding_sequences = gather_cds_files(cds_file_paths)

msa_dir = "01_data/02_phylogenetic_analysis/fasta/04_sco_protein_alignments_v2/"

#prot2transcript_df = pd.read_csv("01_data/02_phylogenetic_analysis/prot2transcript.tsv", 
#   sep = "\t", header = None)
#prot2transcript_df = prot2transcript_df.set_index(0)
#prot2transcript = prot2transcript_df.to_dict(orient = "index")

try:
   os.mkdir("01_data/02_phylogenetic_analysis/fasta/05_cds_for_alignments_v2/")
except:
   None

msa_dir_content = os.listdir(msa_dir)
for i in tqdm(range(len(msa_dir_content)), desc = "Processing MSAs:"):
   msa = msa_dir_content[i]
   get_cds_for_msa(msa_dir + msa, coding_sequences)

