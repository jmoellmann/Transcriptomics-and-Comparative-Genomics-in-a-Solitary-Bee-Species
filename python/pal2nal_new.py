import re
import os
from Bio import Seq
from Bio import SeqRecord
from Bio import SeqIO
import multiprocessing as mp
from tqdm import tqdm


os.chdir("/media/jannik/Samsung_T5/masters_project/01_main_analysis/")

def convert_pal2nal(protein_alignment_fasta, cds_fasta, out_dir):
   orthogroup_accession = os.path.split(protein_alignment_path)[-1]
   aligned_prot_seqs = [seq_rec for seq_rec in SeqIO.parse(protein_alignment_fasta, "fasta")]
   unaligned_cds_seqs = [seq_rec for seq_rec in SeqIO.parse(cds_fasta, "fasta")]
   out_seqs = []
   for prot_seq, cds_seq in zip(aligned_prot_seqs, unaligned_cds_seqs):
      species_name = cds_seq.name.split("|")[0][:-4]
      cds_seq_name = re.search("[XNY]P_[0-9]+\\.[0-9]", cds_seq.name.split("|")[1])[0]
      if len(cds_seq.seq) % 3 != 0:
         print("Length of CDS Sequence", cds_seq_name, "in species", species_name, 
         "not a multiple of three! Length:", len(cds_seq.seq))
         continue
      if cds_seq_name != prot_seq.name:
         print(cds_seq_name, "and", prot_seq.name, "do not match in species", species_name)
         continue
      out_str = Seq.Seq("")
      j = 0
      for i in range(len(prot_seq.seq)):
         if(prot_seq.seq[i] == "-"):
            out_str = out_str + "---"
         else:
            out_str = out_str + cds_seq.seq[j:j+3]
            j += 3
      out_seq = SeqRecord.SeqRecord(out_str, id = species_name + "_" + prot_seq.name, 
      description = "")
      out_seqs.append(out_seq)
   SeqIO.write(out_seqs, out_dir + orthogroup_accession, "fasta")
   
protein_alignment_dir = "01_data/02_phylogenetic_analysis/fasta/04_sco_protein_alignments/"
unaligned_cds_dir = "01_data/02_phylogenetic_analysis/fasta/05_cds_for_alignments/"
out_dir = "01_data/02_phylogenetic_analysis/fasta/06_cds_alignments_test/"

try:
   os.mkdir(out_dir)
except:
   None

protein_alignment_file_paths = [protein_alignment_dir + i for i in os.listdir(protein_alignment_dir)]
unaligned_cds_file_paths = [unaligned_cds_dir + i for i in os.listdir(unaligned_cds_dir)]


for i in tqdm(range(len(protein_alignment_file_paths)), desc = "Processing MSAs"):
   protein_alignment_path = protein_alignment_file_paths[i]
   cds_path = unaligned_cds_file_paths[i]
   convert_pal2nal(protein_alignment_path, cds_path, out_dir)

