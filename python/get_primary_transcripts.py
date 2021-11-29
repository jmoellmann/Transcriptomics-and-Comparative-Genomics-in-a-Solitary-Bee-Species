import os
import HTSeq
import pandas as pd

def get_primary_transcripts(in_fasta, out_fasta, gene2prot):
	gene2prot = pd.read_csv(gene2prot, names = ["GeneID", "ProteinID"])
	gene2prot['GeneID'] = gene2prot['GeneID'].str.replace('GeneID:', '')
	gene2prot['ProteinID'] = gene2prot['ProteinID'].str.replace('Genbank:', '')
	
	sequences = dict( (seq.name, seq) for seq in HTSeq.FastaReader(in_fasta) )
	out_fasta_file = open(out_fasta, "w")
	
	for geneID in gene2prot['GeneID'].unique():
		proteins = gene2prot.loc[gene2prot['GeneID'] == geneID, 'ProteinID'].tolist()
		longest_protein = proteins[0]
		longest_protein_len = len(sequences[proteins[0]])
		for protein in proteins[1:]:
			if len(sequences[protein]) > longest_protein_len:
				longest_protein = protein
				longest_protein_len = len(sequences[protein])
		sequences[longest_protein].write_to_fasta_file(out_fasta_file)
	out_fasta_file.close()
	
	out_sequences = [seq.name for seq in HTSeq.FastaReader(out_fasta)]

	
	print(os.path.split(in_fasta)[1], ":", len(sequences), "proteins read in;", len(gene2prot['ProteinID'].unique()), "protein accessions in gff;", len(out_sequences), "primary transcripts written out;", len(gene2prot['GeneID'].unique()), "gene accessions in gff")

os.chdir("/media/jannik/Samsung_T5/masters_project/ApisTranscriptomics/data/reference/phylogenetic_analysis")
try:
	os.mkdir("fasta/03_primary_transcripts_from_cds")
except:
	None

get_primary_transcripts("fasta/02_proteins_from_cds/Acyrthosiphon_pisum.faa", "fasta/03_primary_transcripts_from_cds/Acyrthosiphon_pisum.faa", "gene2prot/Acyrthosiphon_pisum.txt")
get_primary_transcripts("fasta/02_proteins_from_cds/Aedes_aegypti.faa", "fasta/03_primary_transcripts_from_cds/Aedes_aegypti.faa", "gene2prot/Aedes_aegypti.txt")
get_primary_transcripts("fasta/02_proteins_from_cds/Anopheles_gambiae.faa", "fasta/03_primary_transcripts_from_cds/Anopheles_gambiae.faa", "gene2prot/Anopheles_gambiae.txt")
get_primary_transcripts("fasta/02_proteins_from_cds/Apis_mellifera.faa", "fasta/03_primary_transcripts_from_cds/Apis_mellifera.faa", "gene2prot/Apis_mellifera.txt")
get_primary_transcripts("fasta/02_proteins_from_cds/Bombus_terrestris.faa", "fasta/03_primary_transcripts_from_cds/Bombus_terrestris.faa", "gene2prot/Bombus_terrestris.txt")
get_primary_transcripts("fasta/02_proteins_from_cds/Bombyx_mori.faa", "fasta/03_primary_transcripts_from_cds/Bombyx_mori.faa", "gene2prot/Bombyx_mori.txt")
get_primary_transcripts("fasta/02_proteins_from_cds/Ceratina_calcarata.faa", "fasta/03_primary_transcripts_from_cds/Ceratina_calcarata.faa", "gene2prot/Ceratina_calcarata.txt")
get_primary_transcripts("fasta/02_proteins_from_cds/Drosophila_melanogaster.faa", "fasta/03_primary_transcripts_from_cds/Drosophila_melanogaster.faa", "gene2prot/Drosophila_melanogaster.txt")
get_primary_transcripts("fasta/02_proteins_from_cds/Dufourea_novaeangliae.faa", "fasta/03_primary_transcripts_from_cds/Dufourea_novaeangliae.faa", "gene2prot/Dufourea_novaeangliae.txt")
get_primary_transcripts("fasta/02_proteins_from_cds/Eufriesea_mexicana.faa", "fasta/03_primary_transcripts_from_cds/Eufriesea_mexicana.faa", "gene2prot/Eufriesea_mexicana.txt")
get_primary_transcripts("fasta/02_proteins_from_cds/Habropoda_laboriosa.faa", "fasta/03_primary_transcripts_from_cds/Habropoda_laboriosa.faa", "gene2prot/Habropoda_laboriosa.txt")
get_primary_transcripts("fasta/02_proteins_from_cds/Megachile_rotundata.faa", "fasta/03_primary_transcripts_from_cds/Megachile_rotundata.faa", "gene2prot/Megachile_rotundata.txt")
get_primary_transcripts("fasta/02_proteins_from_cds/Megalopta_genalis.faa", "fasta/03_primary_transcripts_from_cds/Megalopta_genalis.faa", "gene2prot/Megalopta_genalis.txt")
get_primary_transcripts("fasta/02_proteins_from_cds/Nasonia_vitripennis.faa", "fasta/03_primary_transcripts_from_cds/Nasonia_vitripennis.faa", "gene2prot/Nasonia_vitripennis.txt")
get_primary_transcripts("fasta/02_proteins_from_cds/Nomia_melanderi.faa", "fasta/03_primary_transcripts_from_cds/Nomia_melanderi.faa", "gene2prot/Nomia_melanderi.txt")
get_primary_transcripts("fasta/02_proteins_from_cds/Osmia_bicornis.faa", "fasta/03_primary_transcripts_from_cds/Osmia_bicornis.faa", "gene2prot/Osmia_bicornis.txt")
get_primary_transcripts("fasta/02_proteins_from_cds/Osmia_lignaria.faa", "fasta/03_primary_transcripts_from_cds/Osmia_lignaria.faa", "gene2prot/Osmia_lignaria.txt")
get_primary_transcripts("fasta/02_proteins_from_cds/Polistes_dominula.faa", "fasta/03_primary_transcripts_from_cds/Polistes_dominula.faa", "gene2prot/Polistes_dominula.txt")
get_primary_transcripts("fasta/02_proteins_from_cds/Solenopsis_invicta.faa", "fasta/03_primary_transcripts_from_cds/Solenopsis_invicta.faa", "gene2prot/Solenopsis_invicta.txt")
get_primary_transcripts("fasta/02_proteins_from_cds/Tribolium_castaneum.faa", "fasta/03_primary_transcripts_from_cds/Tribolium_castaneum.faa", "gene2prot/Tribolium_castaneum.txt")
get_primary_transcripts("fasta/02_proteins_from_cds/Vespa_mandarinia.faa", "fasta/03_primary_transcripts_from_cds/Vespa_mandarinia.faa", "gene2prot/Vespa_mandarinia.txt")
