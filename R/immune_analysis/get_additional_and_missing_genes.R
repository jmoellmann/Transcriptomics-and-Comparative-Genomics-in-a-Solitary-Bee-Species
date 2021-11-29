library(tidyverse)
library(biomaRt)

project_dir <- "/media/jannik/Samsung_T5/masters_project/ApisTranscriptomics/"
setwd(project_dir)

vectorized_grepl <- Vectorize(grepl, vectorize.args = "pattern")

bt_mart <- biomaRt::useEnsemblGenomes(biomart = "metazoa_mart", dataset = "bterrestris_eg_gene")
dm_mart <- biomaRt::useEnsemblGenomes(biomart = "metazoa_mart", dataset = "dmelanogaster_eg_gene")

bt_imm_genes <- read.table(
   "data/reference/immune_genes/putative_b_terrestris_immune_genes_ensembl.txt")[["V1"]]
dm_imm_genes <- read.table(
   "data/reference/immune_genes/Dm_immune_genes_flybase.txt")[["V1"]]

putative_ob_immune_genes <- read.table("results/PRJNA285788/immune_gene_identification/immune_genes/ob_merged_imm_genes.txt")[[1]]

orthologues <- read_tsv("results/PRJNA285788/immune_gene_identification/orthofinder_insects/Phylogenetic_Hierarchical_Orthogroups/N0.tsv") %>% 
   dplyr::select(-HOG, -"Gene Tree Parent Clade") %>% 
   dplyr::rename(Orthogroup = OG) %>% 
   pivot_longer(cols = Acyrthosiphon_pisum:Vespa_mandarinia, names_to = "Species") %>% 
   dplyr::rename(Proteins = value) %>% mutate(Species = str_replace(Species, "_", " ")) %>% 
   distinct() %>% group_by(Orthogroup, Species) %>% 
   summarise(Orthogroup, Species, Proteins = paste(Proteins, collapse = ", ")) %>% 
   distinct()

bt_gene2bt_protein <- biomaRt::getBM(c("ensembl_gene_id", "ensembl_peptide_id"), 
                                     mart = bt_mart) %>% 
   as_tibble() %>% 
   dplyr::filter(ensembl_peptide_id != "") %>% 
   dplyr::rename("Bt_ProteinID" = "ensembl_peptide_id", "Bt_GeneID" = "ensembl_gene_id")

bt_gene2gene_description <- read_tsv("data/reference/annotation/bt_gene2gene_description.tsv", 
                                     col_names = c("Bt_GeneID", "Bt_GeneDescription"))

dm_gene2dm_protein <- biomaRt::getBM(c("ensembl_gene_id", "refseq_peptide"), 
                                     mart = dm_mart) %>%
   as_tibble() %>% 
   dplyr::filter(refseq_peptide != "") %>% 
   dplyr::rename("Dm_ProteinID" = "refseq_peptide", "Dm_GeneID" = "ensembl_gene_id")

# dm_gene2gene_description <- read_tsv("data/reference/annotation/dm_gene2gene_description.tsv",
#                                      col_names = c("Dm_GeneID_non_flybase", "Dm_GeneDescription"))
# 
# dm_flybase_gene2dm_gene <- read_tsv("data/reference/annotation/dm_fb_gene2dm_gene.tsv",
#                                     col_names = c("Dm_GeneID_non_flybase", "Dm_GeneID"))

dm_flybase_gene2gene_description <- read_tsv("data/reference/annotation/dm_fb_gene2gene_description.tsv", 
                                             col_names = c("Dm_GeneID", "Dm_GeneDescription"))

ob_gene2protein <- read_tsv("data/reference/annotation/ob_gene2protein.tsv", 
                               col_names = c("Ob_GeneID", "Ob_ProteinID"))

ob_gene2gene_description <- read_tsv("data/reference/annotation/ob_gene2gene_description.tsv",
                                     col_names = c("Ob_GeneID", "Ob_GeneDescription"))

ob_protein2prot_length <- read_csv("data/reference/annotation/ob_primary_transcript_lengths.csv") %>% 
   dplyr::select(ProteinID, Length) %>% dplyr::rename(ProteinLength = Length) %>% 
   mutate(ProteinLength = as.character(ProteinLength))

ob_gene2chr_start_end <- read_tsv("data/reference/annotation/ob_gene2chr_start_end.tsv", 
                                  col_names = c("geneID", "chromosome", "strand", "start", "end")) %>% 
   mutate(chromosome = ifelse(str_detect(chromosome, "CAJRAL"), "unplaced", chromosome)) %>% 
   mutate_all(as.character)

ob_gene2interproscan <- as_tibble(read.csv(
   "results/PRJNA285788/immune_gene_identification/interproscan_results/Osmia_bicornis.faa.tsv", 
   header = FALSE, sep = "\t", 
   col.names =  c("ProteinAccession", "SequenceMD5Digest", "SequenceLength", 
                  "Analysis", "SignatureAccession", "SignatureDescription", 
                  "StartLocation", "StopLocation", "Score", "Status", "Date", 
                  "InterProAccession", "InterProDescription"))) %>% 
   dplyr::filter(Score != "-") %>% mutate(Score = as.numeric(Score)) %>% 
   dplyr::filter(Score < 0.05, InterProAccession != "-") %>% 
   inner_join(ob_gene2protein, by = c("ProteinAccession" = "Ob_ProteinID")) %>% 
   dplyr::select(Ob_GeneID, InterProDescription) %>% group_by(Ob_GeneID) %>% 
   summarise(Ob_GeneID, InterProDomains = paste(unique(InterProDescription), collapse = ", ")) %>% 
   distinct()

ol_gene2protein <- read_tsv("data/reference/annotation/ol_gene2protein.tsv", 
                               col_names = c("Ol_GeneID", "Ol_ProteinID"))

ol_gene2gene_description <- read_tsv("data/reference/annotation/ol_gene2gene_description.tsv",
                                     col_names = c("Ol_GeneID", "Ol_GeneDescription"))

ol_protein2prot_length <- read_csv("data/reference/annotation/ol_primary_transcript_lengths.csv") %>% 
   dplyr::select(ProteinID, Length) %>% dplyr::rename(ProteinLength = Length) %>% 
   mutate(ProteinLength = as.character(ProteinLength))

ol_gene2chr_start_end <- read_tsv("data/reference/annotation/ol_gene2chr_start_end.tsv", 
                                  col_names = c("geneID", "chromosome", "strand", "start", "end")) %>% 
   mutate_all(as.character)

dm_imm_prots <- as.vector(unlist(dm_gene2dm_protein %>% 
                                    dplyr::filter(Dm_GeneID %in% dm_imm_genes, Dm_ProteinID != "") %>% 
                                    dplyr::select(Dm_ProteinID) %>% distinct()))

bt_imm_prots <- as.vector(unlist(bt_gene2bt_protein %>% 
                                    dplyr::filter(Bt_GeneID %in% bt_imm_genes, Bt_ProteinID != "") %>% 
                                    dplyr::select(Bt_ProteinID) %>% distinct()))
dm_imm_prot_ogs <- read_tsv(
   "results/PRJNA285788/immune_gene_identification/orthofinder_insects/Phylogenetic_Hierarchical_Orthogroups/N0.tsv") %>% 
   dplyr::select(OG, Drosophila_melanogaster) %>% 
   filter(apply(vectorized_grepl(dm_imm_prots, Drosophila_melanogaster), MARGIN = 1, 
                function(x) any(x == TRUE))) %>% pull(OG) %>% unique()

bt_imm_prot_ogs <- read_tsv(
   "results/PRJNA285788/immune_gene_identification/orthofinder_insects/Phylogenetic_Hierarchical_Orthogroups/N0.tsv") %>% 
   dplyr::select(OG, Bombus_terrestris) %>% 
   filter(apply(vectorized_grepl(bt_imm_prots, Bombus_terrestris), MARGIN = 1, 
                function(x) any(x == TRUE))) %>% pull(OG) %>% unique()

imm_prot_ogs <- union(dm_imm_prot_ogs, bt_imm_prot_ogs)

# Immune genes missing in Osmia ---------------------------------------------------------------

dm_ogs_missing_in_osmia <- orthologues %>% filter(Orthogroup %in% dm_imm_prot_ogs, 
                       Species %in% c("Osmia bicornis", "Osmia lignaria"), Proteins == "NA") %>% 
   group_by(Orthogroup) %>% filter(n() > 1) %>% pull(Orthogroup) %>% unique()

bt_ogs_missing_in_osmia <- orthologues %>% filter(Orthogroup %in% bt_imm_prot_ogs, 
                       Species %in% c("Osmia bicornis", "Osmia lignaria"), Proteins == "NA") %>% 
   group_by(Orthogroup) %>% filter(n() > 1) %>% pull(Orthogroup) %>% unique()

dm_prots_missing_in_osmia <- orthologues %>% 
   filter(Orthogroup %in% dm_ogs_missing_in_osmia, Species == "Drosophila melanogaster", 
          Proteins != "NA") %>% pull(Proteins) %>% paste(collapse = ", ") %>% str_split(", ") %>% 
   unlist() %>% unique()

bt_prots_missing_in_osmia <- orthologues %>% 
   filter(Orthogroup %in% bt_ogs_missing_in_osmia, Species == "Bombus terrestris", 
          Proteins != "NA") %>% pull(Proteins) %>% paste(collapse = ", ") %>% str_split(", ") %>% 
   unlist() %>% unique()

dm_genes_missing_in_osmia <- tibble(Dm_ProteinID = dm_prots_missing_in_osmia) %>% 
   left_join(dm_gene2dm_protein) %>% pull(Dm_GeneID)
dm_genes_missing_in_osmia <- dm_genes_missing_in_osmia[!is.na(dm_genes_missing_in_osmia)]

dm_genes_missing_in_osmia_tbl <- dm_flybase_gene2gene_description %>% 
   filter(Dm_GeneID %in% dm_genes_missing_in_osmia) %>% 
   rename(geneID = Dm_GeneID, geneDescription = Dm_GeneDescription) %>% 
   mutate(origin = "D.melanogaster")

bt_genes_missing_in_osmia <- tibble(Bt_ProteinID = bt_prots_missing_in_osmia) %>% 
   left_join(bt_gene2bt_protein) %>% pull(Bt_GeneID)
bt_genes_missing_in_osmia <- bt_genes_missing_in_osmia[!is.na(bt_genes_missing_in_osmia)]

bt_genes_missing_in_osmia_tbl <- bt_gene2gene_description %>% 
   filter(Bt_GeneID %in% bt_genes_missing_in_osmia) %>% 
   rename(geneID = Bt_GeneID, geneDescription = Bt_GeneDescription) %>% 
   mutate(origin = "B.terrestris")

genes_missing_in_osmia_tbl <- rbind(dm_genes_missing_in_osmia_tbl, bt_genes_missing_in_osmia_tbl) %>% 
   dplyr::select(geneID, origin, geneDescription)

write_tsv(genes_missing_in_osmia_tbl, "results/PRJNA285788/immune_gene_identification/immune_genes_missing_in_osmia.tsv")

# Additional genes in Osmia (not immune-specific) ---------------------------------------------

osmia_shared_ogs <- orthologues %>% filter(Species %in% c("Osmia lignaria", "Osmia bicornis")) %>% 
   filter(Proteins != "NA") %>% group_by(Orthogroup) %>% filter(n() > 1) %>% pull(Orthogroup) %>% 
   unique()

osmia_unique_ogs <- orthologues %>% filter(!(Species %in% c("Osmia bicornis", "Osmia lignaria"))) %>% 
   group_by(Orthogroup) %>% filter(all(Proteins == "NA")) %>% ungroup() %>% 
   filter(Orthogroup %in% osmia_shared_ogs) %>% pull(Orthogroup) %>% unique()

osmia_unique_genes <- orthologues %>% filter(Orthogroup %in% osmia_unique_ogs, Species %in% c("Osmia bicornis", "Osmia lignaria")) %>% 
   separate_rows(Proteins, sep = ", ") %>% arrange(Species, Orthogroup, Proteins) %>% 
   dplyr::select(Species, Orthogroup, Proteins) %>% 
   left_join(ob_gene2protein, by = c("Proteins" = "Ob_ProteinID")) %>% 
   left_join(ol_gene2protein, by = c("Proteins" = "Ol_ProteinID")) %>% 
   unite("GeneID", Ob_GeneID:Ol_GeneID, sep = "") %>% mutate(GeneID = str_remove(GeneID, "NA")) %>% 
   left_join(ob_gene2gene_description, by = c("GeneID" = "Ob_GeneID")) %>% 
   left_join(ol_gene2gene_description, by = c("GeneID" = "Ol_GeneID")) %>% 
   unite("GeneDescription", Ob_GeneDescription:Ol_GeneDescription, sep = "") %>% 
   mutate(GeneDescription = str_remove(GeneDescription, "NA")) %>% 
   left_join(ob_protein2prot_length, by = c("Proteins" = "ProteinID")) %>% 
   left_join(ol_protein2prot_length, by = c("Proteins" = "ProteinID")) %>% 
   unite("ProteinLength", ProteinLength.x:ProteinLength.y, sep = "") %>% 
   mutate(ProteinLength = as.numeric(str_remove(ProteinLength, "NA"))) %>% 
   left_join(ob_gene2chr_start_end, by = c("GeneID" = "geneID")) %>% 
   left_join(ol_gene2chr_start_end, by = c("GeneID" = "geneID")) %>% 
   dplyr::select(Species:ProteinLength, chromosome.x, chromosome.y, strand.x, strand.y, start.x, 
                 start.y, end.x, end.y) %>% 
   unite("Chromosome/Contig", chromosome.x:chromosome.y, sep = "") %>% 
   mutate(`Chromosome/Contig` = str_remove(`Chromosome/Contig`, "NA")) %>% 
   unite("Strand", strand.x:strand.y, sep = "") %>% 
   mutate(Strand = str_remove(Strand, "NA")) %>% 
   unite(Start, start.x:start.y, sep = "") %>% 
   mutate(Start = as.numeric(str_remove(Start, "NA"))) %>% 
   unite(End, end.x:end.y, sep = "") %>% 
   mutate(End = as.numeric(str_remove(End, "NA"))) %>% 
   left_join(ob_gene2interproscan, by = c("GeneID" = "Ob_GeneID")) %>% 
   group_by(Species, Orthogroup) %>% mutate(NOrthologues = n()) %>% 
   na_if("NA") %>% dplyr::select(-Proteins)

write_tsv(osmia_unique_genes, "results/PRJNA285788/immune_gene_identification/osmia_specific_genes.tsv")
