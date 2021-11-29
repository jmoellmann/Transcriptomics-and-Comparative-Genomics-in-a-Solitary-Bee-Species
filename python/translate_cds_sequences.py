import os
import re
import HTSeq
from Bio import Seq
from tqdm import tqdm

os.chdir("/media/jannik/Samsung_T5/masters_project/ApisTranscriptomics/")

def translate_cds_to_protein(cds_path, out_dir):
   species_name = os.path.split(cds_path)[-1][:-4]
   out_file = open(out_dir + species_name + ".faa", "w")
   cds = [(seq.name, seq) for seq in HTSeq.FastaReader(cds_path)]
   for i in tqdm(range(len(cds)), desc = "Processing " + species_name + "...: "):
      sequence_name = cds[i][0]
      seq = cds[i][1]
      biopython_seq = Seq.Seq(str(seq))
      protein_seq = biopython_seq.translate(stop_symbol = "", to_stop = False)
      prot_id = re.search("[XNY]P_[0-9]+\\.[0-9]", sequence_name)[0]
      prot_descr = re.search("\\[protein=.*?\\]", seq.descr)[0].split("=")[1][:-1]
      out_seq = HTSeq.Sequence(str(protein_seq).encode("UTF-8"))
      out_seq.name = prot_id + " " + prot_descr + " [" + species_name + "]"
      out_seq.write_to_fasta_file(out_file)
   out_file.close()

cds_dir = "data/reference/phylogenetic_analysis/fasta/01_cds_from_refseq/"
cds_file_paths = [cds_dir + i for i in os.listdir(cds_dir)]

out_dir = "data/reference/phylogenetic_analysis/fasta/02_proteins_from_cds/"

try:
   os.mkdir(out_dir)
except:
   None

for cds_path in cds_file_paths:
   species_name = os.path.split(cds_path)[-1][:-4]
   if not os.path.exists(out_dir + species_name):
      translate_cds_to_protein(cds_path, out_dir)

