import os
from Bio import pairwise2
from Bio import Seq
import HTSeq
import re
from tqdm import tqdm

os.chdir("/media/jannik/Samsung_T5/masters_project/ApisTranscriptomics/")

def compare_cds_to_protein(protein_path, cds_dict, out_dir):
   species_name = os.path.split(protein_path)[-1][:-4]
   out_file = open(out_dir + species_name + ".txt", "w")
   prot = [(seq.name, seq) for seq in HTSeq.FastaReader(protein_path)]
   for i in tqdm(range(len(prot)), desc = "Processing " + species_name + "...: "):
      sequence_id = prot[i][0]
      seq = prot[i][1]
      biopython_seq = Seq.Seq(str(cds_dict[sequence_id]))
      cds_seq = biopython_seq.translate(stop_symbol = "", to_stop = False)
      prot_seq = Seq.Seq(str(seq))
      if len(prot_seq) != len(cds_seq):
         out_file.write("1>---------------------" + "\n")
         out_file.write(sequence_id + " in " + species_name + " with difference in prot. length / transl. cds length: " + str(len(prot_seq)) + " / " + str(len(cds_seq)) + ":\n")
         out_file.write(str(prot_seq) + "\n" + str(cds_seq) + "\n\n\n")
      elif pairwise2.align.globalxx(cds_seq, prot_seq, score_only = True, one_alignment_only=True) != len(prot_seq):
         out_file.write("2>---------------------" + "\n")
         out_file.write(sequence_id + " in " + species_name + " with non-perfect alignment between prot. and transl. cds:\n")
         out_file.write(str(prot_seq) + "\n" + str(cds_seq) + "\n\n\n")
   out_file.close()
         
def write_cds_to_dict(cds_file_paths):
   coding_sequences = dict()
   for cds_file_path in cds_file_paths:
      species_name = os.path.split(cds_file_path)[-1][:-3]
      temp_dict_v1 = dict( (seq.name, seq) for seq in HTSeq.FastaReader(cds_file_path) )
      temp_dict_v2 = dict()
      for key, value in temp_dict_v1.items():
         temp_key = re.search("[XNY]P_[0-9]+\\.[0-9]", key)[0]
         temp_seq = value
         temp_seq.name = species_name + "_" + temp_seq.name[:-2] + "_" + temp_seq.name[:-1]
         temp_dict_v2[temp_key] = temp_seq
      coding_sequences.update(temp_dict_v2)
   return coding_sequences


cds_dir = "data/reference/phylogenetic_analysis/fasta/cds_from_genomic/"
cds_file_paths = [cds_dir + i for i in os.listdir(cds_dir)]
cds_dict = write_cds_to_dict(cds_file_paths)

species_list = [i[:-3] for i in os.listdir(cds_dir)]

proteins_dir = "data/reference/phylogenetic_analysis/fasta/proteins/"
proteins_file_paths = [proteins_dir + i for i in os.listdir(proteins_dir) if i[:-4] in species_list]

out_dir = "data/reference/phylogenetic_analysis/fasta/cds_protein_comparison/"

try:
   os.mkdir(out_dir)
except:
   None

for protein_path in proteins_file_paths:
   species_name = os.path.split(protein_path)[-1][:-4]
   if not os.path.exists(out_dir + species_name):
      compare_cds_to_protein(protein_path, cds_dict, out_dir)

