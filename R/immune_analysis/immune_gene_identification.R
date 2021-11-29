library(ggplot2)
library(dplyr)
library(stringr)
library(readr)
library(tidyr)
library(biomaRt)
library(eulerr)
library(scales)
library(ggpubr)

get_type <- function(value, bt_imm_prots, dm_imm_prots){
  if(value %in% bt_imm_prots & value %in% dm_imm_prots){
    return("Dm/Bt shared")
  } else if(value %in% bt_imm_prots){
    return("Bt unique")
  } else if(value %in% dm_imm_prots){
    return("Dm unique")
  } else {
    return(NA)
  }
}

project_dir <- "/media/jannik/Samsung_T5/masters_project/ApisTranscriptomics/"
setwd(project_dir)

insect_orthologues_dir <- 
  "results/PRJNA285788/immune_gene_identification/orthofinder_insects/Phylogenetic_Hierarchical_Orthogroups"
bee_orthologues_dir <- 
  "results/PRJNA285788/immune_gene_identification/orthofinder_bees/Phylogenetic_Hierarchical_Orthogroups"

dir.create("results/PRJNA285788/immune_gene_identification/plots/seq_identity", 
           showWarnings = FALSE, recursive = TRUE)
dir.create("results/PRJNA285788/immune_gene_identification/plots/prot_lengths_and_domains", 
           showWarnings = FALSE, recursive = TRUE)
dir.create("results/PRJNA285788/immune_gene_identification/plots/overlap", 
           showWarnings = FALSE, recursive = TRUE)
dir.create("results/PRJNA285788/immune_gene_identification/funct_enrichment", 
           showWarnings = FALSE, recursive = TRUE)
# Load bt immune genes and dm immune peptides -------------------------------------------------

bt_imm_genes <- read.table(
  "data/reference/immune_genes/putative_b_terrestris_immune_genes_ensembl.txt")[["V1"]]
dm_imm_genes <- read.table(
  "data/reference/immune_genes/Dm_immune_genes_flybase.txt")[["V1"]]
am_imm_genes <- paste0("LOC", read.table(
  "data/reference/immune_genes/putative_a_mellifera_immune_genes_2017_ST14.txt")[["V1"]])

# Connect to ensembl metazoa ------------------------------------------------------------------

#bi_mart <- biomaRt::useEnsemblGenomes(biomart = "metazoa_mart", dataset = "bimpatiens_eg_gene")
am_mart <- biomaRt::useEnsemblGenomes(biomart = "metazoa_mart", dataset = "amellifera_eg_gene")
bt_mart <- biomaRt::useEnsemblGenomes(biomart = "metazoa_mart", dataset = "bterrestris_eg_gene")
dm_mart <- biomaRt::useEnsemblGenomes(biomart = "metazoa_mart", dataset = "dmelanogaster_eg_gene")

# Orthofinder protein orthologue tables -------------------------------------------------------

orthofinder_res <- readr::read_tsv(paste0(insect_orthologues_dir, "/N0.tsv")) 
orthofinder_res_bees <- readr::read_tsv(paste0(bee_orthologues_dir, "/N0.tsv")) 

dm_protein2ob_protein_orthofinder <- orthofinder_res %>% 
  dplyr::select(OG, Osmia_bicornis, Drosophila_melanogaster) %>%
  tidyr::separate_rows(Osmia_bicornis, sep = ", ") %>% 
  tidyr::separate_rows(Drosophila_melanogaster, sep = ", ") %>% 
  dplyr::rename("Ob_ProteinID" = "Osmia_bicornis", "Dm_ProteinID" = "Drosophila_melanogaster") %>% 
  dplyr::select("Ob_ProteinID", "Dm_ProteinID") %>% drop_na()

dm_protein2bt_protein_orthofinder <- orthofinder_res %>% 
  dplyr::select(OG, Bombus_terrestris, Drosophila_melanogaster) %>%
  tidyr::separate_rows(Drosophila_melanogaster, sep = ", ") %>% 
  tidyr::separate_rows(Bombus_terrestris, sep = ", ") %>% 
  dplyr::rename("Dm_ProteinID" = "Drosophila_melanogaster", "Bt_ProteinID" = "Bombus_terrestris") %>% 
  dplyr::select("Dm_ProteinID", "Bt_ProteinID") %>% drop_na()

ob_protein2bt_protein_orthofinder <- orthofinder_res_bees %>% 
  dplyr::select(OG, Bombus_terrestris, Osmia_bicornis) %>%
  tidyr::separate_rows(Osmia_bicornis, sep = ", ") %>% 
  tidyr::separate_rows(Bombus_terrestris, sep = ", ") %>% 
  dplyr::rename("Ob_ProteinID" = "Osmia_bicornis", "Bt_ProteinID" = "Bombus_terrestris") %>% 
  dplyr::select("Ob_ProteinID", "Bt_ProteinID")

dm_protein2am_protein_orthofinder <- orthofinder_res %>% 
  dplyr::select(OG, Apis_mellifera, Drosophila_melanogaster) %>%
  tidyr::separate_rows(Drosophila_melanogaster, sep = ", ") %>%
  tidyr::separate_rows(Apis_mellifera, sep = ", ") %>%
  dplyr::rename("Dm_ProteinID" = "Drosophila_melanogaster", "Am_ProteinID" = "Apis_mellifera") %>%
  dplyr::select("Dm_ProteinID", "Am_ProteinID")

ob_protein2am_protein_orthofinder <- orthofinder_res %>% 
  dplyr::select(OG, Apis_mellifera, Osmia_bicornis) %>%
  tidyr::separate_rows(Osmia_bicornis, sep = ", ") %>%
  tidyr::separate_rows(Apis_mellifera, sep = ", ") %>%
  dplyr::rename("Ob_ProteinID" = "Osmia_bicornis", "Am_ProteinID" = "Apis_mellifera") %>%
  dplyr::select("Ob_ProteinID", "Am_ProteinID")

# Gene to Protein tables ----------------------------------------------------------------------

ob_gene2ob_protein <- read_tsv("data/reference/annotation/ob_gene2protein.tsv", 
                               col_names = c("Ob_GeneID", "Ob_ProteinID"))

bt_gene2bt_protein <- biomaRt::getBM(c("ensembl_gene_id", "ensembl_peptide_id"), 
                                     mart = bt_mart) %>% 
  as_tibble() %>% 
  dplyr::filter(ensembl_peptide_id != "") %>% 
  dplyr::rename("Bt_ProteinID" = "ensembl_peptide_id", "Bt_GeneID" = "ensembl_gene_id")

dm_gene2dm_protein <- biomaRt::getBM(c("ensembl_gene_id", "refseq_peptide"), 
                                     mart = dm_mart) %>%
  as_tibble() %>% 
  dplyr::filter(refseq_peptide != "") %>% 
  dplyr::rename("Dm_ProteinID" = "refseq_peptide", "Dm_GeneID" = "ensembl_gene_id")

dm_gene2dm_protein_flybase <- biomaRt::getBM(c("flybase_gene_id", "refseq_peptide"), 
                                     mart = dm_mart) %>%
  as_tibble() %>% 
  dplyr::filter(refseq_peptide != "") %>% 
  dplyr::rename("Dm_ProteinID" = "refseq_peptide", "Dm_GeneID" = "flybase_gene_id")

am_gene2am_protein <- as_tibble(biomaRt::getBM(c("ensembl_gene_id", "ensembl_peptide_id"),
                                               mart = am_mart)) %>%
  dplyr::filter(ensembl_peptide_id != "") %>%
  dplyr::rename("Am_ProteinID" = "ensembl_peptide_id", "Am_GeneID" = "ensembl_gene_id")

am_gene2am_protein_from_gtf <- read_tsv("data/reference/annotation/am_gene2protein.tsv",
                                        col_names = c("Am_GeneID", "Am_ProteinID"))

# bi_gene2bi_protein <- biomaRt::getBM(c("ensembl_gene_id", "refseq_peptide"), 
#                                                   mart = bi_mart) %>% 
#   as_tibble() %>% 
#   dplyr::filter(refseq_peptide != "") %>% 
#   dplyr::rename("Bi_ProteinID" = "refseq_peptide", "Bi_GeneID" = "ensembl_gene_id")
# 

ob_gene2am_gene <- ob_protein2am_protein_orthofinder %>% 
  mutate(Am_ProteinID = str_remove(Am_ProteinID, "\\.[0-9]")) %>% left_join(am_gene2am_protein) %>% 
  left_join(ob_gene2ob_protein) %>% dplyr::select(Am_GeneID, Ob_GeneID) %>% drop_na()

write_csv(ob_gene2am_gene, "results/PRJNA285788/immune_gene_identification/ob_gene2am_gene.csv")

bt_gene2am_gene <- biomaRt::getBM(attributes = c("ensembl_gene_id", "amellifera_eg_homolog_ensembl_gene"), 
               mart = bt_mart) %>% as_tibble() %>% 
  dplyr::rename("Bt_GeneID" = "ensembl_gene_id", 
                "Am_GeneID" = "amellifera_eg_homolog_ensembl_gene") %>% 
  dplyr::filter(Am_GeneID != "", Bt_GeneID != "")

bt_am_imm_genes <- bt_gene2am_gene %>% dplyr::filter(Am_GeneID %in% am_imm_genes) %>% 
  dplyr::select(Bt_GeneID) %>% unlist() %>% as.vector()

# Infer Osmia immune genes --------------------------------------------------------------------

dm_imm_prots <- as.vector(unlist(dm_gene2dm_protein %>% 
                                   dplyr::filter(Dm_GeneID %in% dm_imm_genes, Dm_ProteinID != "") %>% 
                                   dplyr::select(Dm_ProteinID) %>% distinct()))

bt_imm_prots <- as.vector(unlist(bt_gene2bt_protein %>% 
                                   dplyr::filter(Bt_GeneID %in% bt_imm_genes, Bt_ProteinID != "") %>% 
                                   dplyr::select(Bt_ProteinID) %>% distinct()))

am_imm_prots <- as.vector(unlist(am_gene2am_protein %>% 
                                   dplyr::filter(Am_GeneID %in% am_imm_genes, Am_ProteinID != "") %>% 
                                   dplyr::select(Am_ProteinID) %>% distinct()))

dm_bt_imm_prots <- unique(as.vector(unlist(
  dm_protein2bt_protein_orthofinder %>% dplyr::filter(Bt_ProteinID %in% bt_imm_prots) %>% 
    dplyr::select(Dm_ProteinID))))
ob_bt_imm_prots <- unique(as.vector(unlist(
  ob_protein2bt_protein_orthofinder %>% dplyr::filter(Bt_ProteinID %in% bt_imm_prots) %>% 
    dplyr::select(Ob_ProteinID))))
ob_dm_imm_prots <- unique(as.vector(unlist(
  dm_protein2ob_protein_orthofinder %>% dplyr::filter(Dm_ProteinID %in% dm_imm_prots) %>% 
    dplyr::select(Ob_ProteinID))))
bt_dm_imm_prots <- unique(as.vector(unlist(
  dm_protein2bt_protein_orthofinder %>% dplyr::filter(Dm_ProteinID %in% dm_imm_prots) %>% 
    dplyr::select(Bt_ProteinID))))

ob_bt_imm_genes <- unique(as.vector(unlist(
  left_join(ob_gene2ob_protein, left_join(ob_protein2bt_protein_orthofinder, 
                                          bt_gene2bt_protein)) %>% 
    filter(!is.na(Bt_GeneID), Bt_GeneID %in% bt_imm_genes) %>% 
    dplyr::select(Ob_GeneID))))

ob_dm_imm_genes <- unique(as.vector(unlist(
  left_join(ob_gene2ob_protein, dm_protein2ob_protein_orthofinder) %>% 
    filter(!is.na(Dm_ProteinID), Dm_ProteinID %in% dm_imm_prots) %>%
    dplyr::select(Ob_GeneID))))


ob_am_imm_genes <- ob_gene2am_gene %>% dplyr::filter(Am_GeneID %in% am_imm_genes) %>% 
  pull(Ob_GeneID)

bt_dm_imm_genes <- unique(as.vector(unlist(
  left_join(bt_gene2bt_protein, dm_protein2bt_protein_orthofinder) %>% 
    filter(!is.na(Dm_ProteinID), Dm_ProteinID %in% dm_imm_prots) %>% 
    dplyr::select(Bt_GeneID))))

dm_bt_imm_genes <- unique(as.vector(unlist(
  left_join(dm_gene2dm_protein, dm_protein2bt_protein_orthofinder) %>% 
    filter(!is.na(Bt_ProteinID), Bt_ProteinID %in% bt_imm_prots) %>% 
    dplyr::select(Dm_GeneID))))

dm_am_imm_genes <- unique(as.vector(unlist(
  left_join(dm_gene2dm_protein, dm_protein2am_protein_orthofinder) %>% 
    mutate(Am_ProteinID = str_remove(Am_ProteinID, pattern = "\\.[0-9]$")) %>% 
    filter(!is.na(Am_ProteinID), Am_ProteinID %in% am_imm_prots) %>% 
    dplyr::select(Dm_GeneID))))

dm_merged_imm_genes <- unique(c(dm_imm_genes, dm_bt_imm_genes))
bt_merged_imm_genes <- unique(c(bt_imm_genes, bt_dm_imm_genes))
ob_merged_imm_genes <- unique(c(ob_dm_imm_genes, ob_bt_imm_genes))

dm_merged_imm_prots <- unique(c(dm_imm_prots, dm_bt_imm_prots))
ob_merged_imm_prots <- unique(c(ob_dm_imm_prots, ob_bt_imm_prots))
bt_merged_imm_prots <- unique(c(bt_dm_imm_prots, bt_imm_prots))

dir.create("results/PRJNA285788/immune_gene_identification/immune_genes", recursive = TRUE, 
           showWarnings = FALSE)

write.table(bt_imm_genes, "results/PRJNA285788/immune_gene_identification/immune_genes/bt_imm_genes.txt", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(dm_imm_genes, "results/PRJNA285788/immune_gene_identification/immune_genes/dm_imm_genes.txt", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(am_imm_genes, "results/PRJNA285788/immune_gene_identification/immune_genes/am_imm_genes.txt", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(dm_am_imm_genes, "results/PRJNA285788/immune_gene_identification/immune_genes/dm_am_imm_genes.txt", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(ob_bt_imm_genes, "results/PRJNA285788/immune_gene_identification/immune_genes/ob_bt_imm_genes.txt", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(ob_dm_imm_genes, "results/PRJNA285788/immune_gene_identification/immune_genes/ob_dm_imm_genes.txt", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(ob_am_imm_genes, "results/PRJNA285788/immune_gene_identification/immune_genes/ob_am_imm_genes.txt", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(bt_dm_imm_genes, "results/PRJNA285788/immune_gene_identification/immune_genes/bt_dm_imm_genes.txt", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(dm_bt_imm_genes, "results/PRJNA285788/immune_gene_identification/immune_genes/dm_bt_imm_genes.txt", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(dm_merged_imm_genes, "results/PRJNA285788/immune_gene_identification/immune_genes/dm_merged_imm_genes.txt", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(bt_merged_imm_genes, "results/PRJNA285788/immune_gene_identification/immune_genes/bt_merged_imm_genes.txt", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(ob_merged_imm_genes, "results/PRJNA285788/immune_gene_identification/immune_genes/ob_merged_imm_genes.txt", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(ob_merged_imm_prots, "results/PRJNA285788/immune_gene_identification/immune_genes/ob_merged_imm_prots.txt", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(bt_merged_imm_prots, "results/PRJNA285788/immune_gene_identification/immune_genes/bt_merged_imm_prots.txt", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(dm_merged_imm_prots, "results/PRJNA285788/immune_gene_identification/immune_genes/dm_merged_imm_prots.txt", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)


# Protein Lengths -----------------------------------------------------------------------------

dm_protein_lengths <- read_csv("data/reference/annotation/dm_primary_transcript_lengths.csv") %>% 
  dplyr::select(ProteinID, Length) %>% dplyr::rename("Dm_ProteinID" = "ProteinID")
dm_protein_lengths_imm_genes <- dm_protein_lengths %>% 
  dplyr::filter(Dm_ProteinID %in% dm_merged_imm_prots)

ob_protein_lengths <- read_csv("data/reference/annotation/ob_primary_transcript_lengths.csv") %>% 
  dplyr::select(ProteinID, Length) %>% dplyr::rename("Ob_ProteinID" = "ProteinID")
ob_protein_lengths_imm_genes <- ob_protein_lengths %>% 
  dplyr::filter(Ob_ProteinID %in% ob_merged_imm_prots)

bt_protein_lengths <- read_csv("data/reference/annotation/bt_primary_transcript_lengths.csv") %>% 
  dplyr::select(ProteinID, Length) %>% dplyr::rename("Bt_ProteinID" = "ProteinID")

am_protein_lengths <- read_csv("data/reference/annotation/am_primary_transcript_lengths.csv") %>% 
  dplyr::select(ProteinID, Length) %>% dplyr::rename("Am_ProteinID" = "ProteinID")

ob_dm_bt_am_prot_lengths <- tibble(Species = c(rep("O.bicornis", length(ob_protein_lengths$Length)), 
                                               rep("D.melanogaster", length(dm_protein_lengths$Length)),
                                               rep("B.terrestris", length(bt_protein_lengths$Length)),
                                               rep("A.mellifera", length(am_protein_lengths$Length))),
                                   ProteinLength = c(ob_protein_lengths$Length, 
                                                     dm_protein_lengths$Length,
                                                     bt_protein_lengths$Length,
                                                     am_protein_lengths$Length))

ob_dm_bt_am_prot_lengths_imm_genes <- tibble(Species = c(rep("O.bicornis", length(ob_protein_lengths_imm_genes$Length)), 
                                                 rep("D.melanogaster", length(dm_protein_lengths_imm_genes$Length))),
                                     ProteinLength = c(ob_protein_lengths_imm_genes$Length, 
                                                       dm_protein_lengths_imm_genes$Length))

(prot_length_distributions_plot <- ggplot(ob_dm_bt_am_prot_lengths) + 
    geom_density(aes(x = ProteinLength, fill = Species), alpha = .5) + 
    xlim(0, 1800) + 
    ggtitle("Distributions of protein lengths of \nproteins in O.bicornis and D.melanogaster"))

(prot_length_imm_genes_distributions_plot <- ggplot(ob_dm_bt_am_prot_lengths_imm_genes) + 
    geom_density(aes(x = ProteinLength, fill = Species), alpha = .5) + 
    xlim(0, 1800) + 
    ggtitle("Distributions of protein lengths of \nimmune proteins in O.bicornis and D.melanogaster"))

ggsave("results/PRJNA285788/immune_gene_identification/plots/overlap/dm_ob_bt_am_prot_len_distr.png", 
       prot_length_distributions_plot)
ggsave("results/PRJNA285788/immune_gene_identification/plots/overlap/dm_ob_imm_prot_len_distr.png", 
       prot_length_imm_genes_distributions_plot)

#t.test(ob_protein_lengths$Length, dm_protein_lengths$Length)
t.test(ob_protein_lengths_imm_genes$Length, dm_protein_lengths_imm_genes$Length)


# Get protein orthologue pair Sequence identity ---------------------------------------------

sequenceIDs_insects <- read_tsv("results/PRJNA285788/immune_gene_identification/orthofinder_insects/WorkingDirectory/SequenceIDs.tsv", 
                        col_names = c("Accession", "EnsemblProteinID", "Description")) %>% 
  mutate(Accession = str_remove(Accession, ":"))

sequenceIDs_bees <- read_tsv("results/PRJNA285788/immune_gene_identification/orthofinder_bees/WorkingDirectory/SequenceIDs.tsv", 
                             col_names = c("Accession", "EnsemblProteinID", "Description")) %>% 
  mutate(Accession = str_remove(Accession, ":"))


dm2ob_stats <- read_tsv("results/PRJNA285788/immune_gene_identification/orthofinder_insects/WorkingDirectory/Blast7_15.txt", 
                        col_names = c("QueryAccession", "TargetAccession", "SequenceIdentity", "AlignmentLength", 
                                      "V5", "V6", "V7", "V8", "V9", "V10", "EValue", "V12")) %>% 
  dplyr::select(-(V5:V10), -V12) %>% 
  left_join(sequenceIDs_insects, by = c("QueryAccession" = "Accession")) %>% 
  left_join(sequenceIDs_insects, by = c("TargetAccession" = "Accession"), suffix = c("_dm", "_ob")) %>% 
  dplyr::select(EnsemblProteinID_dm, EnsemblProteinID_ob, SequenceIdentity, AlignmentLength, EValue) %>% 
  dplyr::rename(Dm_ProteinID = EnsemblProteinID_dm, Ob_ProteinID = EnsemblProteinID_ob) %>% 
  left_join(ob_protein_lengths) %>% mutate(SequenceSimilarity = AlignmentLength / Length * SequenceIdentity) %>%
  right_join(dm_protein2ob_protein_orthofinder) %>% 
  dplyr::arrange(Dm_ProteinID, Ob_ProteinID) %>% 
  mutate(Type = ifelse(Ob_ProteinID %in% ob_merged_imm_prots, 
                                "immune-protein", "non-immune-protein")) %>% drop_na()  %>% 
  dplyr::mutate(Type2 = sapply(Dm_ProteinID, get_type, dm_bt_imm_prots, dm_imm_prots)) %>% 
  dplyr::mutate(Type3 = ifelse(Type2 %in% c("Dm unique", "Dm/Bt shared"), "Dm", "Non-Dm"))

bt2ob_stats <- read_tsv("results/PRJNA285788/immune_gene_identification/orthofinder_bees/WorkingDirectory/Blast1_9.txt", 
                        col_names = c("QueryAccession", "TargetAccession", "SequenceIdentity", "AlignmentLength", 
                                      "V5", "V6", "V7", "V8", "V9", "V10", "EValue", "V12")) %>% 
  dplyr::select(-(V5:V10), -V12) %>% 
  left_join(sequenceIDs_bees, by = c("QueryAccession" = "Accession")) %>% 
  left_join(sequenceIDs_bees, by = c("TargetAccession" = "Accession"), suffix = c("_bt", "_ob")) %>% 
  dplyr::select(EnsemblProteinID_bt, EnsemblProteinID_ob, SequenceIdentity, AlignmentLength, EValue) %>% 
  dplyr::rename(Bt_ProteinID = EnsemblProteinID_bt, Ob_ProteinID = EnsemblProteinID_ob) %>% 
  left_join(ob_protein_lengths) %>% mutate(SequenceSimilarity = AlignmentLength / Length * SequenceIdentity) %>%
  right_join(ob_protein2bt_protein_orthofinder) %>% 
  dplyr::arrange(Bt_ProteinID, Ob_ProteinID) %>% 
  mutate(Type = ifelse(Ob_ProteinID %in% ob_merged_imm_prots, 
                                "immune-protein", "non-immune-protein")) %>% drop_na() %>% 
  dplyr::mutate(Type2 = sapply(Bt_ProteinID, get_type, bt_imm_prots, bt_dm_imm_prots)) %>% 
  dplyr::mutate(Type3 = ifelse(Type2 %in% c("Bt unique", "Dm/Bt shared"), "Bt", "Non-Bt"))



dm2ob_stats_imm_prots <- dm2ob_stats %>% dplyr::filter(Ob_ProteinID %in% ob_merged_imm_prots)
dm2ob_stats_dm_imm_prots <- dm2ob_stats %>% dplyr::filter(Dm_ProteinID %in% dm_imm_prots)
dm2ob_stats_bt_imm_prots <- dm2ob_stats %>% dplyr::filter(Dm_ProteinID %in% dm_bt_imm_prots)
dm2ob_stats_dm_only_imm_prots <- dm2ob_stats %>% dplyr::filter(Dm_ProteinID %in% dm_imm_prots, 
                                                               !(Dm_ProteinID %in% dm_bt_imm_prots))
dm2ob_stats_bt_only_imm_prots <- dm2ob_stats %>% dplyr::filter(Dm_ProteinID %in% dm_bt_imm_prots, 
                                                               !(Dm_ProteinID %in% dm_imm_prots))
dm2ob_stats_dm_bt_shared_imm_prots <- dm2ob_stats %>% dplyr::filter(Dm_ProteinID %in% dm_imm_prots, 
                                                               Dm_ProteinID %in% dm_bt_imm_prots)
bt2ob_stats_imm_prots <- bt2ob_stats %>% dplyr::filter(Ob_ProteinID %in% ob_merged_imm_prots) %>% 
  dplyr::mutate(Type2 = sapply(Bt_ProteinID, get_type, bt_imm_prots, bt_dm_imm_prots))
bt2ob_stats_bt_imm_prots <- bt2ob_stats %>% dplyr::filter(Bt_ProteinID %in% bt_imm_prots)

dm2ob_stats_bt_dm_imm_prots <- tibble(Type = c(rep("B.terrestris unique", 
                                                               length(dm2ob_stats_bt_only_imm_prots$SequenceIdentity)),
                                                    rep("D.melanogaster unique",
                                                        length(dm2ob_stats_dm_only_imm_prots$SequenceIdentity)),
                                                    rep("B.t./D.m. shared",
                                                    length(dm2ob_stats_dm_bt_shared_imm_prots$SequenceIdentity))),
                                                    SequenceIdentity = c(dm2ob_stats_bt_only_imm_prots$SequenceIdentity,
                                                                           dm2ob_stats_dm_only_imm_prots$SequenceIdentity,
                                                                           dm2ob_stats_dm_bt_shared_imm_prots$SequenceIdentity))
bt_dm_imm_prot_seq_ident <- tibble(Type = c(rep("B.terrestris derived", 
                                                   length((bt2ob_stats %>% dplyr::filter(Type3 == "Bt"))$SequenceIdentity)),
                                               rep("D.melanogaster derived",
                                                   length((dm2ob_stats %>% dplyr::filter(Type3 == "Dm"))$SequenceIdentity))),
                                      SequenceIdentity = c((bt2ob_stats %>% dplyr::filter(Type3 == "Bt"))$SequenceIdentity,
                                                             (dm2ob_stats %>% dplyr::filter(Type3 == "Dm"))$SequenceIdentity))

# ob2dm_stats <- read_tsv("results/PRJNA285788/immune_gene_identification/orthofinder_insects/WorkingDirectory/Blast15_7.txt", 
#                         col_names = c("QueryAccession", "TargetAccession", "SequenceIdentity", 
#                                       "AlignmentLength", 
#                                       "V5", "V6", "V7", "V8", "V9", "V10", "EValue", "V12")) %>% 
#   dplyr::select(-(V5:V10), -V12) %>% 
#   left_join(sequenceIDs, by = c("QueryAccession" = "Accession")) %>% 
#   left_join(sequenceIDs, by = c("TargetAccession" = "Accession"), suffix = c("_ob", "_dm")) %>% 
#   dplyr::select(EnsemblProteinID_ob, EnsemblProteinID_dm, SequenceIdentity, AlignmentLength, EValue) %>% 
#   rename(Dm_ProteinID = EnsemblProteinID_dm, Ob_ProteinID = EnsemblProteinID_ob) %>% 
#   left_join(dm_protein_lengths) %>% mutate(SequenceIdentity = AlignmentLength / Length * SequenceIdentity) %>% 
#   dplyr::arrange(Dm_ProteinID, Ob_ProteinID) %>% 
#   drop_na()

(dm2ob_seq_ident_distr_plot <- ggplot(dm2ob_stats) +
    geom_density(aes(x = SequenceIdentity, fill = Type), alpha = .5) + 
    ggtitle("Sequence identity between O.bicornis and D.melanogaster protein orthologues"))
# (dm2ob_bt_dm_imm_prot_seq_identity_distr_plot <- ggplot(dm2ob_stats_bt_dm_imm_prots) +
#   geom_density(aes(x = SequenceIdentity, fill = Type), alpha = .5) + 
#     ggtitle("Distributions of Sequence identity between O.bicornis and D.melanogaster"))
(bt2ob_seq_ident_distr_plot <- ggplot(bt2ob_stats) +
    geom_density(aes(x = SequenceIdentity, fill = Type), alpha = .5) + 
    ggtitle("Sequence identity between O.bicornis and B.terrestris protein orthologues"))

(bt_dm_imm_gene_seq_ident_distr_plot <- ggplot(bt_dm_imm_prot_seq_ident) +
    geom_density(aes(x = SequenceIdentity, fill = Type), alpha = .65) + 
    ggtitle("Sequence identity of O.bicornis immune proteins derived via B.terrestris/D.melanogaster"))

(unique_imm_gene_bt_seq_ident_distr_plot <- ggplot(drop_na(bt2ob_stats)) +
    geom_density(aes(x = SequenceIdentity, fill = Type2), alpha = 0.65) + 
    xlab("Sequence identity [%]") + ylab("Density") + 
    scale_fill_discrete(name = "Immune\nprotein origin", labels = c("B.t. only", "D.m. only", 
                                                                    "B.t. + D.m.")))

(unique_imm_gene_dm_seq_ident_distr_plot <- ggplot(drop_na(dm2ob_stats)) +
    geom_density(aes(x = SequenceIdentity, fill = Type2), alpha = .5) + 
    xlab("Sequence identity [%]") + ylab("Density") + 
    scale_fill_discrete(name = "Immune\nprotein origin", 
                        labels = c("B.t. only", "D.m. only", "B.t. + D.m.")))


(bt2ob_seq_ident_vs_length_plot <- ggplot(bt2ob_stats, aes(x = Length, y = SequenceIdentity,
                                                      colour = Type)) +
    geom_point(alpha = 0.1) + geom_smooth(method = "glm") + xlim(0,2000) + 
    ggtitle("Protein length vs Sequence identity for O.bicornis protein orthologues (via B.terrestris)"))
(dm2ob_seq_ident_vs_length_plot <- ggplot(dm2ob_stats, aes(x = Length, y = SequenceIdentity,
                                                      colour = Type)) +
    geom_point(alpha = 0.1) + geom_smooth(method = "glm") + xlim(0,2000) + 
    ggtitle("Protein length vs Sequence identity for O.bicornis protein orthologues (via D.melanogaster)"))

cor.test(bt2ob_stats[bt2ob_stats$Type == "immune-protein", ]$Length,
         bt2ob_stats[bt2ob_stats$Type == "immune-protein", ]$SequenceIdentity)
cor.test(bt2ob_stats[bt2ob_stats$Type == "non-immune-protein", ]$Length,
         bt2ob_stats[bt2ob_stats$Type == "non-immune-protein", ]$SequenceIdentity)
cor.test(dm2ob_stats[dm2ob_stats$Type == "immune-protein", ]$Length,
         dm2ob_stats[dm2ob_stats$Type == "immune-protein", ]$SequenceIdentity)
cor.test(dm2ob_stats[dm2ob_stats$Type == "non-immune-protein", ]$Length,
         dm2ob_stats[dm2ob_stats$Type == "non-immune-protein", ]$SequenceIdentity)

ggsave("results/PRJNA285788/immune_gene_identification/plots/seq_identity/dm2ob_distr.png", 
       dm2ob_seq_ident_distr_plot, width = 25, height = 15, units = "cm")
ggsave("results/PRJNA285788/immune_gene_identification/plots/seq_identity/bt2ob_distr.png", 
       bt2ob_seq_ident_distr_plot, width = 25, height = 15, units = "cm")
ggsave("results/PRJNA285788/immune_gene_identification/plots/seq_identity/bt_dm_imm_gene_distr.png", 
       bt_dm_imm_gene_seq_ident_distr_plot, width = 25, height = 15, units = "cm")
ggsave("results/PRJNA285788/immune_gene_identification/plots/seq_identity/bt2ob_seq_ident_vs_length.png", 
       bt2ob_seq_ident_vs_length_plot, width = 25, height = 15, units = "cm")
ggsave("results/PRJNA285788/immune_gene_identification/plots/seq_identity/dm2ob_seq_ident_vs_length.png", 
       dm2ob_seq_ident_vs_length_plot, width = 25, height = 15, units = "cm")
ggsave("results/PRJNA285788/immune_gene_identification/plots/seq_identity/unique_imm_gene_bt_seq_ident_distr.png", 
       unique_imm_gene_bt_seq_ident_distr_plot, width = 25, height = 15, units = "cm")
ggsave("results/PRJNA285788/immune_gene_identification/plots/seq_identity/unique_imm_gene_dm_seq_ident_distr.png", 
       unique_imm_gene_dm_seq_ident_distr_plot, width = 25, height = 15, units = "cm")



t.test(dm2ob_stats$SequenceIdentity, dm2ob_stats_imm_prots$SequenceIdentity)
t.test(dm2ob_stats_dm_only_imm_prots$SequenceIdentity, dm2ob_stats_bt_imm_prots$SequenceIdentity)
t.test(bt2ob_stats$SequenceIdentity, bt2ob_stats_imm_prots$SequenceIdentity)

# Ensembl protein orthologue tables -----------------------------------------------------------

dm_refseq_protein2dm_ensembl_protein <- tibble::as_tibble(
  biomaRt::getBM(c("ensembl_peptide_id", "refseq_peptide"), mart = dm_mart)) %>% 
  dplyr::rename("Dm_ProteinID_Ensembl" = "ensembl_peptide_id", 
                "Dm_ProteinID_RefSeq" = "refseq_peptide") 

dm_protein2bt_protein_ensembl <- tibble::as_tibble(
  biomaRt::getBM(c("ensembl_peptide_id", "bterrestris_eg_homolog_ensembl_peptide"), mart = dm_mart)) %>% 
  dplyr::rename("Dm_ProteinID_Ensembl" = "ensembl_peptide_id", 
                "Bt_ProteinID" = "bterrestris_eg_homolog_ensembl_peptide") %>% 
  dplyr::left_join(dm_refseq_protein2dm_ensembl_protein) %>% 
  dplyr::rename("Dm_ProteinID" = "Dm_ProteinID_RefSeq") %>% 
  dplyr::filter(Dm_ProteinID != "", !is.na(Dm_ProteinID), 
                str_detect(Dm_ProteinID, "NP_[0-9]+\\.[0-9]"),
                str_detect(Bt_ProteinID, "XP_[0-9]+\\.[0-9]")) %>% 
  dplyr::select(Dm_ProteinID, Bt_ProteinID) %>% arrange(Dm_ProteinID) %>% dplyr::distinct()

dm_protein2am_protein_ensembl <- tibble::as_tibble(
  biomaRt::getBM(c("ensembl_peptide_id", "amellifera_eg_homolog_ensembl_peptide"), mart = dm_mart)) %>% 
  dplyr::rename("Dm_ProteinID_Ensembl" = "ensembl_peptide_id", 
                "Am_ProteinID" = "amellifera_eg_homolog_ensembl_peptide") %>% 
  dplyr::left_join(dm_refseq_protein2dm_ensembl_protein) %>% 
  dplyr::rename("Dm_ProteinID" = "Dm_ProteinID_RefSeq") %>% 
  dplyr::filter(Dm_ProteinID != "", !is.na(Dm_ProteinID), 
                str_detect(Dm_ProteinID, "NP_[0-9]+\\.[0-9]"),
                str_detect(Am_ProteinID, "XP_[0-9]+")) %>% 
  dplyr::select(Dm_ProteinID, Am_ProteinID) %>% arrange(Dm_ProteinID) %>% dplyr::distinct()



dm_ensembl_gene2dm_flybase_annotation <- biomaRt::getBM(
  attributes = c("ensembl_gene_id", "flybase_annotation_id"), mart = dm_mart) %>% as_tibble() %>% 
  dplyr::filter(str_detect(flybase_annotation_id, "CG")) %>% 
  dplyr::rename("Dm_GeneID" = ensembl_gene_id, "Dm_FlybaseID" = flybase_annotation_id)

dm_imm_genes_2013 <- read.table("data/reference/annotation/dm_immune_gene_list_2013.txt", 
                                header = FALSE, col.names = c("Dm_FlybaseID")) %>% 
  left_join(dm_ensembl_gene2dm_flybase_annotation) %>% drop_na() %>% dplyr::select(Dm_GeneID) %>% 
  unlist() %>% as.vector()

(overlap_imm_genes_w2013_plot <- plot(euler(list(B.terrestris = dm_bt_imm_genes, 
                                          D.melanogaster = dm_imm_genes,
                                      D.melanogaster_2013 = dm_imm_genes_2013)), 
                                quantities = TRUE, 
                                main = "Overlap of D.melanogaster and B.terrestris immune genes"))
(overlap_imm_genes_wAm_plot <- plot(euler(list(B.terrestris = dm_bt_imm_genes, 
                                           A.mellifera = dm_am_imm_genes,
                                                 D.melanogaster = dm_imm_genes)), 
                                      quantities = TRUE, 
                                      main = "Overlap of D.melanogaster and B.terrestris/A.mellifera immune genes"))
(overlap_imm_genes_plot <- plot(euler(list(B.terrestris = ob_bt_imm_genes, 
                                               D.melanogaster = ob_dm_imm_genes)), 
                                    quantities = TRUE, 
                                    main = "Overlap of D.melanogaster and B.terrestris immune genes"))

(overlap_imm_genes_plot_no_title <- plot(euler(list(`B. terrestris` = ob_bt_imm_genes, 
                                           `D. melanogaster` = ob_dm_imm_genes)), 
                                quantities = TRUE))

ggsave("results/PRJNA285788/immune_gene_identification/plots/overlap//imm_gene_overlap_w2013.png", 
       overlap_imm_genes_w2013_plot, width = 25, height = 15, units = "cm")
ggsave("results/PRJNA285788/immune_gene_identification/plots/overlap//imm_gene_overlap.png", 
       overlap_imm_genes_plot, width = 25, height = 15, units = "cm")
ggsave("results/PRJNA285788/immune_gene_identification/plots/overlap/imm_gene_wAm_overlap.png", 
       overlap_imm_genes_wAm_plot, width = 25, height = 15, units = "cm")

# Gene to Gene Description --------------------------------------------------------------------

bt_gene2gene_description <- read_tsv("data/reference/annotation/bt_gene2gene_description.tsv", 
                                 col_names = c("Bt_GeneID", "Bt_GeneDescription")) %>% 
  arrange(Bt_GeneID)

# Overlap between ensembl Dm/Bt orthologues and OrthoFinder Dm/Bt orthologues -----------------

dm_gene2bt_gene_orthofinder <- dplyr::left_join(
  dplyr::right_join(dm_protein2bt_protein_orthofinder, dm_gene2dm_protein), bt_gene2bt_protein) %>% 
  dplyr::filter(Dm_GeneID != "", !is.na(Dm_GeneID), Bt_GeneID != "", !is.na(Bt_GeneID)) %>% 
  dplyr::select(Dm_GeneID, Bt_GeneID) %>% dplyr::distinct() 

dm_gene2bt_gene_ensembl <- dplyr::left_join(
  dplyr::right_join(dm_protein2bt_protein_ensembl, dm_gene2dm_protein), bt_gene2bt_protein) %>% 
  dplyr::filter(Dm_GeneID != "", !is.na(Dm_GeneID), Bt_GeneID != "", !is.na(Bt_GeneID)) %>% 
  dplyr::select(Dm_GeneID, Bt_GeneID) %>% dplyr::distinct()

orthl_pairs_prot <- list()
orthl_pairs_prot[["OrthoFinder"]] <- as.vector(unlist(tidyr::unite(
  dm_protein2bt_protein_orthofinder, Orthologue_Relation, Dm_ProteinID, Bt_ProteinID) %>% 
                                              arrange(Orthologue_Relation)))
orthl_pairs_prot[["Ensembl"]] <- as.vector(unlist(tidyr::unite(
  dm_protein2bt_protein_ensembl, Orthologue_Relation, Dm_ProteinID, Bt_ProteinID) %>% 
                                         arrange(Orthologue_Relation)))

orthl_pairs_gene <- list()
orthl_pairs_gene[["OrthoFinder"]] <- as.vector(unlist(tidyr::unite(
  dm_gene2bt_gene_orthofinder, Orthologue_Relation, Dm_GeneID, Bt_GeneID) %>% 
                                                            arrange(Orthologue_Relation)))
orthl_pairs_gene[["Ensembl"]] <- as.vector(unlist(tidyr::unite(
  dm_gene2bt_gene_ensembl, Orthologue_Relation, Dm_GeneID, Bt_GeneID) %>% 
                                                        arrange(Orthologue_Relation)))

orthl_genes <- list()
orthl_genes[["Drosophila"]][["Ensembl"]] <- unique(dm_gene2bt_gene_ensembl$Dm_GeneID)
orthl_genes[["Drosophila"]][["OrthoFinder"]] <- unique(dm_gene2bt_gene_orthofinder$Dm_GeneID)
orthl_genes[["Drosophila"]][["Dmelanogaster_imm_genes"]] <- as.vector(unlist(dm_gene2dm_protein %>% 
  dplyr::filter(Dm_ProteinID %in% dm_imm_prots) %>% 
  dplyr::select(Dm_GeneID) %>% distinct()))
orthl_genes[["Bombus"]][["Ensembl"]] <- unique(dm_gene2bt_gene_ensembl$Bt_GeneID)
orthl_genes[["Bombus"]][["OrthoFinder"]] <- unique(dm_gene2bt_gene_orthofinder$Bt_GeneID)
orthl_genes[["Bombus"]][["Bterrestris_imm_genes"]] <- bt_imm_genes

prot_lvl_comparison_plot <- (plot(euler(orthl_pairs_prot), quantities = TRUE, 
                            main = "B.terrestris/D.melanogaster orthologue protein pairs inferred via Ensembl and OrthoFinder"))
gene_lvl_comparison_plot <- plot(euler(orthl_pairs_gene), quantities = TRUE,
                            main = "B.terrestris/D.melanogaster orthologue genes pairs inferred via Ensembl and OrthoFinder")
orthl_genes_drosophila_plot <- plot(euler(orthl_genes[["Drosophila"]]), quantities = TRUE,
     main = "D.melanogaster genes for which a B. terrestris orthologue exists in Ensembl/OrthoFinder")
orthl_genes_bombus_plot <- plot(euler(orthl_genes[["Bombus"]]), quantities = TRUE,
     main = "B.terrestris genes for which a D. melanogaster orthologue exists in Ensembl/OrthoFinder")

ggsave("results/PRJNA285788/immune_gene_identification/plots/overlap/prot_orthl_pair_overlap.png", 
       prot_lvl_comparison_plot, width = 32, height = 18, units = "cm")
ggsave("results/PRJNA285788/immune_gene_identification/plots/overlap/gene_orthl_pair_overlap.png", 
       gene_lvl_comparison_plot, width = 32, height = 18, units = "cm")
ggsave("results/PRJNA285788/immune_gene_identification/plots/overlap/genes_drosophila_plot.png", 
       orthl_genes_drosophila_plot, width = 32, height = 18, units = "cm")
ggsave("results/PRJNA285788/immune_gene_identification/plots/overlap/genes_bombus_plot.png", 
       orthl_genes_bombus_plot, width = 32, height = 18, units = "cm")


# Interpro domain annotation ------------------------------------------------------------------

dm_gene2interpro_ensembl <- as_tibble(biomaRt::getBM(c("interpro", "interpro_description", 
                                               "ensembl_gene_id"), mart = dm_mart)) %>% 
  dplyr::filter(ensembl_gene_id %in% dm_imm_genes, interpro != "-") %>% 
  dplyr::rename("Dm_GeneID" = "ensembl_gene_id", "InterproID" = "interpro", 
                "InterproDescription" = "interpro_description")

bt_gene2interpro_ensembl <- as_tibble(biomaRt::getBM(c("interpro", "interpro_description", 
                                               "ensembl_gene_id"), mart = bt_mart)) %>% 
  dplyr::filter(interpro != "-") %>% 
  dplyr::rename("Bt_GeneID" = "ensembl_gene_id", "InterproID" = "interpro", 
                "InterproDescription" = "interpro_description")
bt_gene2interpro_ensembl_dm_imm_genes <- bt_gene2interpro_ensembl %>% 
  dplyr::filter(Bt_GeneID %in% bt_dm_imm_genes)

bt_gene2interproscan <- as_tibble(read.csv(
  "results/PRJNA285788/immune_gene_identification/interproscan_results/Bombus_terrestris.faa.tsv", 
  header = FALSE, sep = "\t", 
  col.names =  c("ProteinAccession", "SequenceMD5Digest", "SequenceLength", 
                 "Analysis", "SignatureAccession", "SignatureDescription", 
                 "StartLocation", "StopLocation", "Score", "Status", "Date", 
                 "InterProAccession", "InterProDescription"))) %>% 
  dplyr::filter(Score != "-") %>% mutate(Score = as.numeric(Score)) %>% 
  dplyr::filter(Score < 0.05, InterProAccession != "-") %>% 
  inner_join(bt_gene2bt_protein, by = c("ProteinAccession" = "Bt_ProteinID"))
bt_gene2interproscan_dm_imm_genes <- bt_gene2interproscan %>% 
  dplyr::filter(Bt_GeneID %in% bt_dm_imm_genes)

print(paste0("BT: ", length(unique(bt_gene2interproscan$ProteinAccession)) / 13240))

ob_gene2interproscan <- as_tibble(read.csv(
  "results/PRJNA285788/immune_gene_identification/interproscan_results/Osmia_bicornis.faa.tsv", 
  header = FALSE, sep = "\t", 
  col.names =  c("ProteinAccession", "SequenceMD5Digest", "SequenceLength", 
                 "Analysis", "SignatureAccession", "SignatureDescription", 
                 "StartLocation", "StopLocation", "Score", "Status", "Date", 
                 "InterProAccession", "InterProDescription"))) %>% 
  dplyr::filter(Score != "-") %>% mutate(Score = as.numeric(Score)) %>% 
  dplyr::filter(Score < 0.05, InterProAccession != "-") %>% 
  inner_join(ob_gene2ob_protein, by = c("ProteinAccession" = "Ob_ProteinID"))
ob_gene2interproscan_dm_imm_genes <- ob_gene2interproscan %>% 
  dplyr::filter(Ob_GeneID %in% ob_dm_imm_genes)

print(paste0("OB: ", length(unique(ob_gene2interproscan$ProteinAccession)) / 13851))

dm_gene2interproscan <- as_tibble(read.csv(
  "results/PRJNA285788/immune_gene_identification/interproscan_results/Drosophila_melanogaster.faa.tsv", 
  header = FALSE, sep = "\t", 
  col.names =  c("ProteinAccession", "SequenceMD5Digest", "SequenceLength", 
                 "Analysis", "SignatureAccession", "SignatureDescription", 
                 "StartLocation", "StopLocation", "Score", "Status", "Date", 
                 "InterProAccession", "InterProDescription"))) %>% 
  dplyr::filter(Score != "-") %>% mutate(Score = as.numeric(Score)) %>% 
  dplyr::filter(Score < 0.05, InterProAccession != "-") %>% 
  inner_join(dm_gene2dm_protein, by = c("ProteinAccession" = "Dm_ProteinID"))
dm_gene2interproscan_dm_imm_genes <- dm_gene2interproscan %>% 
  dplyr::filter(Dm_GeneID %in% dm_imm_genes)

print(paste0("DM: ", length(unique(dm_gene2interproscan$ProteinAccession)) / 13969))

ensembl_vs_interproscan_dm_imm_genes <- list()
ensembl_vs_interproscan_dm_imm_genes[["InterPro Ensembl"]] <- unique(
  bt_gene2interpro_ensembl_dm_imm_genes$InterproID)
ensembl_vs_interproscan_dm_imm_genes[["InterProScan"]] <- unique(
  bt_gene2interproscan_dm_imm_genes$InterProAccession)

ensembl_vs_interproscan_dm_imm_genes_plot <- plot(eulerr::euler(ensembl_vs_interproscan_dm_imm_genes), 
                                       quantities = TRUE, 
                                       main = "Number of InterPro domains inferred via InterProScan/Ensembl in B.terrestris")

ggsave("results/PRJNA285788/immune_gene_identification/plots/overlap/domain_overlap_interproscan_ensembl_bt.png", 
       ensembl_vs_interproscan_dm_imm_genes_plot, width = 32, height = 18, units = "cm")

interproscan_dm_imm_genes <- list()
interproscan_dm_imm_genes[["B.terrestris"]] <- unique(
  bt_gene2interproscan_dm_imm_genes$InterProAccession)
interproscan_dm_imm_genes[["D.melanogaster"]] <- unique(
  dm_gene2interproscan_dm_imm_genes$InterProAccession)
interproscan_dm_imm_genes[["O.bicornis"]] <- unique(
  ob_gene2interproscan_dm_imm_genes$InterProAccession)

# unique_bt_imm_domains <- unique(bt_gene2interpro_dm_imm_genes[!(bt_gene2interpro_dm_imm_genes$InterproDescription %in% 
#                                          dm_gene2interpro_dm_imm_genes$InterproDescription), ])
# 
# write_csv(unique_bt_imm_domains, 
#           "results/PRJNA285788/immune_gene_identification/plots/overlap/unique_bt_imm_domains.csv", 
#           col_names = TRUE)

interproscan_dm_imm_genes_plot <- plot(eulerr::venn(interproscan_dm_imm_genes), quantities = TRUE,  
                                   fills = c("thistle3", "moccasin", "lightsteelblue"),
  main = "Overlap of interproscan domains between Drosophila immune genes and Bombus and Osmia orthologues")

ggsave("results/PRJNA285788/immune_gene_identification/plots/overlap/overlap_imm_gene_orthl_domains.png", 
       interproscan_dm_imm_genes_plot, width = 32, height = 18, units = "cm")

# Gene to GO tables ---------------------------------------------------------------------------

bt_gene2bt_GO <- biomaRt::getBM(attributes = c("ensembl_gene_id", "go_id"), mart = bt_mart) %>% 
  as_tibble() %>% dplyr::rename("Bt_GeneID" = "ensembl_gene_id", "Bt_GOID" = "go_id") %>% 
  dplyr::filter("Bt_GeneID" != "", !is.na("Bt_GeneID"), "Bt_GOID" != "", !is.na("Bt_GOID")) %>% 
  tidyr::drop_na()

dm_gene2dm_GO <- biomaRt::getBM(attributes = c("ensembl_gene_id", "go_id"), mart = dm_mart) %>% 
  as_tibble() %>% dplyr::rename("Dm_GeneID" = "ensembl_gene_id", "Dm_GOID" = "go_id") %>% 
  dplyr::filter("Dm_GeneID" != "" & !is.na("Dm_GeneID"), "Dm_GOID" != "", !is.na("Dm_GOID")) %>% 
  tidyr::drop_na()

ob_gene2dm_GO <- dplyr::left_join(ob_gene2ob_protein,
                                  dplyr::left_join(dm_protein2ob_protein_orthofinder,
                                                   dplyr::left_join(dm_gene2dm_protein,
                                                                    dm_gene2dm_GO))) %>%
  filter(!is.na(Dm_GOID)) %>% dplyr::select("Ob_GeneID", "Dm_GOID")

ob_gene2bt_GO <- dplyr::left_join(ob_gene2ob_protein,
                                  dplyr::left_join(ob_protein2bt_protein_orthofinder,
                                                   dplyr::left_join(bt_gene2bt_protein,
                                                                    bt_gene2bt_GO))) %>%
  filter(!is.na(Bt_GOID)) %>% dplyr::select("Ob_GeneID", "Bt_GOID")

ob_gene2dm_gene <- dplyr::left_join(ob_gene2ob_protein,
                                    dplyr::left_join(dm_protein2ob_protein_orthofinder,
                                                     dm_gene2dm_protein)) %>%
  filter(!is.na(Dm_GeneID)) %>% dplyr::select("Ob_GeneID", "Dm_GeneID")

ob_gene2bt_gene <- dplyr::left_join(ob_gene2ob_protein,
                                    dplyr::left_join(ob_protein2bt_protein_orthofinder,
                                                     bt_gene2bt_protein)) %>%
  filter(!is.na(Bt_GeneID)) %>% dplyr::select("Ob_GeneID", "Bt_GeneID")

ob_gene2dm_gene_flybase <- dplyr::left_join(ob_gene2ob_protein,
                                            dplyr::left_join(dm_protein2ob_protein_orthofinder,
                                                             dm_gene2dm_protein_flybase)) %>%
  filter(!is.na(Dm_GeneID)) %>% dplyr::select("Ob_GeneID", "Dm_GeneID")

write_csv(ob_gene2dm_GO, "results/PRJNA285788/immune_gene_identification/funct_enrichment/ob_gene2dm_GO.csv")
write_csv(ob_gene2bt_GO, "results/PRJNA285788/immune_gene_identification/funct_enrichment/ob_gene2bt_GO.csv")
write_csv(ob_gene2dm_gene, "results/PRJNA285788/immune_gene_identification/ob_gene2dm_gene.csv")
write_csv(ob_gene2dm_gene_flybase, "results/PRJNA285788/immune_gene_identification/ob_gene2dm_gene_fb.csv")

# Domain type overlap per gene ----------------------------------------------------------------

dm_gene2unique_interpro_domains <- dm_gene2interproscan %>% group_by(Dm_GeneID) %>% 
  summarise(Dm_domains = paste(unique(InterProAccession), collapse = ", "))

ob_gene2unique_interpro_domains <- ob_gene2interproscan %>% group_by(Ob_GeneID) %>% 
  summarise(Ob_domains = paste(unique(InterProAccession), collapse = ", "))

bt_gene2unique_interpro_domains <- bt_gene2interproscan %>% group_by(Bt_GeneID) %>% 
  summarise(Bt_domains = paste(unique(InterProAccession), collapse = ", "))

vectorized_intersect <- Vectorize(intersect, vectorize.args = c("x", "y"))
vectorized_unique <- Vectorize(unique, vectorize.args = "x")

ob_dm_domain_overlap <- ob_gene2dm_gene %>% left_join(dm_gene2unique_interpro_domains) %>% 
  left_join(ob_gene2unique_interpro_domains) %>% drop_na() %>% 
  mutate(Overlap = sapply(vectorized_intersect(str_split(Dm_domains, pattern = ", "), 
                                               str_split(Ob_domains, pattern = ", ")), FUN = length)) %>% 
  mutate(Total = sapply(vectorized_unique(str_split(paste(Dm_domains, Ob_domains, sep = ", "), pattern = ", ")), FUN = length)) %>% 
  mutate(RelativeOverlap = Overlap / Total)

ob_dm_domain_overlap_imm_genes <- ob_dm_domain_overlap %>% dplyr::filter(Ob_GeneID %in% ob_merged_imm_genes)

ob_bt_domain_overlap <- ob_gene2bt_gene %>% left_join(bt_gene2unique_interpro_domains) %>% 
  left_join(ob_gene2unique_interpro_domains) %>% drop_na() %>% 
  mutate(Overlap = sapply(vectorized_intersect(str_split(Bt_domains, pattern = ", "), 
                                               str_split(Ob_domains, pattern = ", ")), FUN = length)) %>% 
  mutate(Total = sapply(vectorized_unique(str_split(paste(Bt_domains, Ob_domains, sep = ", "), pattern = ", ")), FUN = length)) %>% 
  mutate(RelativeOverlap = Overlap / Total)

ob_bt_domain_overlap_imm_genes <- ob_bt_domain_overlap %>% dplyr::filter(Ob_GeneID %in% ob_merged_imm_genes)

# Number of domains per protein ---------------------------------------------------------------

dm_prot_len2bt_prot_len <- left_join(left_join(dm_protein_lengths, dm_protein2bt_protein_orthofinder), 
                                     bt_protein_lengths, 
                                     by = "Bt_ProteinID", suffix = c("_Dm", "_Bt")) %>% 
  dplyr::filter(!is.na(Length_Dm), !is.na(Length_Bt), Dm_ProteinID %in% dm_imm_prots)

dm_prot_len2ob_prot_len <- left_join(left_join(dm_protein_lengths, dm_protein2ob_protein_orthofinder), 
                                     ob_protein_lengths, 
                                     by = "Ob_ProteinID", suffix = c("_Dm", "_Ob")) %>% 
  dplyr::filter(!is.na(Length_Dm), !is.na(Length_Ob), Ob_ProteinID %in% ob_merged_imm_prots)

bt_prot_len2ob_prot_len <- left_join(left_join(bt_protein_lengths, ob_protein2bt_protein_orthofinder), 
                                     ob_protein_lengths, 
                                     by = "Ob_ProteinID", suffix = c("_Bt", "_Ob")) %>% 
  dplyr::filter(!is.na(Length_Bt), !is.na(Length_Ob), 
                Bt_ProteinID %in% dm_protein2bt_protein_orthofinder[
                  dm_protein2bt_protein_orthofinder$Dm_ProteinID %in% dm_imm_prots, ]$Bt_ProteinID)

dm_prot2num_domains <- dm_gene2interproscan %>% 
  group_by(ProteinAccession) %>% 
  summarise(Dm_Domains = list(sort(unique(InterProAccession)))) %>% 
  mutate(Dm_NumDomains = lengths(Dm_Domains)) %>% 
  dplyr::rename(Dm_ProteinID = ProteinAccession)
bt_prot2num_domains <- bt_gene2interproscan %>% 
  group_by(ProteinAccession) %>% 
  summarise(Bt_Domains = list(sort(unique(InterProAccession)))) %>% 
  mutate(Bt_NumDomains = lengths(Bt_Domains)) %>% 
  dplyr::rename(Bt_ProteinID = ProteinAccession)
ob_prot2num_domains <- ob_gene2interproscan %>% 
  group_by(ProteinAccession) %>% 
  summarise(Ob_Domains = list(sort(unique(InterProAccession)))) %>% 
  mutate(Ob_NumDomains = lengths(Ob_Domains)) %>% 
  dplyr::rename(Ob_ProteinID = ProteinAccession)

# dm_ob_bt_num_domains <- tibble(Species = c(rep("D.melanogaster", length(dm_prot2num_domains$Dm_NumDomains)),
#                                            rep("B.terrestris", length(bt_prot2num_domains$Bt_NumDomains)),
#                                            rep("O.bicornis", length(ob_prot2num_domains$Ob_NumDomains))),
#                                NumDomains = c(dm_prot2num_domains$Dm_NumDomains,
#                                               bt_prot2num_domains$Bt_NumDomains,
#                                               ob_prot2num_domains$Ob_NumDomains))
# 
# 
# (prot_length_distributions_plot <- ggplot(dm_ob_bt_num_domains) + 
#     geom_density(aes(x = NumDomains, fill = Species), alpha = .5) + 
#     ggtitle("Distributions of protein lengths of \nproteins in O.bicornis and D.melanogaster"))


dm_stats2bt_stats <- inner_join(left_join(dm_prot2num_domains, dm_protein2bt_protein_orthofinder), 
                                bt_prot2num_domains) %>% 
  mutate(DiffNumDomains = Dm_NumDomains - Bt_NumDomains, 
         NumOverlapDomains = mapply(function(x, y) length(intersect(x, y)), Dm_Domains, Bt_Domains), 
         AllEqualDomains = (NumOverlapDomains == Dm_NumDomains & NumOverlapDomains == Bt_NumDomains)) %>% 
  right_join(dm_prot_len2bt_prot_len) %>% drop_na()

dm_stats2ob_stats <- inner_join(left_join(dm_prot2num_domains, dm_protein2ob_protein_orthofinder), 
                                ob_prot2num_domains) %>% 
  mutate(DiffNumDomains = Dm_NumDomains - Ob_NumDomains, 
         NumOverlapDomains = mapply(function(x, y) length(intersect(x, y)), Dm_Domains, Ob_Domains), 
         AllEqualDomains = (NumOverlapDomains == Dm_NumDomains & NumOverlapDomains == Ob_NumDomains),
         ProteinOrigin = ifelse((Ob_ProteinID %in% ob_bt_imm_prots & 
                                  !(Ob_ProteinID %in% ob_dm_imm_prots)), "B.t. only", "D.m. only")) %>%
  mutate(ProteinOrigin = ifelse((ProteinOrigin == "D.m. only") & (Ob_ProteinID %in% ob_bt_imm_prots),
                                "B.t. + D.m.", ProteinOrigin)) %>%
  mutate(ProteinOrigin = factor(ProteinOrigin, ordered = TRUE,
                                   levels = c( "B.t. only", "D.m. only", "B.t. + D.m."))) %>% 
  right_join(dm_prot_len2ob_prot_len) %>% drop_na() %>% 
  dplyr::filter(Ob_ProteinID != "XP_029047272.1")

bt_stats2ob_stats <- inner_join(left_join(bt_prot2num_domains, ob_protein2bt_protein_orthofinder), 
                                ob_prot2num_domains) %>% 
  mutate(DiffNumDomains = Bt_NumDomains - Ob_NumDomains, 
         NumOverlapDomains = mapply(function(x, y) length(intersect(x, y)), Bt_Domains, Ob_Domains), 
         AllEqualDomains = (NumOverlapDomains == Bt_NumDomains & NumOverlapDomains == Ob_NumDomains)) %>% 
  right_join(bt_prot_len2ob_prot_len) %>% drop_na() 


# Protein Length comparison plots -------------------------------------------------------------

(dm_stats2bt_stats_plot <- ggplot(dm_stats2bt_stats) + 
  geom_point(aes(x = Length_Bt, y = Length_Dm, colour = DiffNumDomains)) + 
  scale_colour_viridis_c(name = "Domain \ncount \ndifference", limits = c(-6, 6)) + 
  ylab("D.melanogaster protein length") + 
  xlab("B.terrestris protein length") + 
  xlim(0, max(c(dm_stats2bt_stats$Length_Dm, 
                dm_stats2bt_stats$Length_Bt) * 1.05)) + 
  ylim(0, max(c(dm_stats2bt_stats$Length_Dm, 
                dm_stats2bt_stats$Length_Bt) * 1.05)))
  
(dm_stats2bt_stats_plot_log10 <- ggplot(dm_stats2bt_stats) + 
  geom_point(aes(x = Length_Bt, y = Length_Dm, colour = DiffNumDomains)) + 
  scale_colour_viridis_c(name = "Domain \ncount \ndifference", limits = c(-6, 6)) + 
  scale_x_log10() + 
  scale_y_log10() + ylab("D.melanogaster protein length (log10-transformed)") + 
  xlab("B.terrestris protein length (log10-transformed)"))

(dm_stats2ob_stats_plot_v1 <- ggplot(dm_stats2ob_stats, aes(x = Length_Ob, y = Length_Dm, 
                                                         colour = DiffNumDomains)) + 
  geom_point(alpha = 0.8, size = 1) + geom_smooth(method = "lm") + stat_cor() + 
  scale_colour_viridis_c(name = "Diff. in\ndomain\ncount\n", limits = c(-6, 6)) +  
  ylab(expression(paste(italic("D. melanogaster"),  " protein length [aa]"))) + 
  xlab(expression(paste(italic("O. bicornis"), " protein length [aa]")))) 

(dm_stats2ob_stats_plot_v2 <- ggplot(dm_stats2ob_stats, aes(x = Length_Ob, y = Length_Dm)) + 
    geom_point(alpha = 0.8, size = 1, aes(colour = ProteinOrigin)) + geom_smooth(method = "lm") + stat_cor() + 
    scale_colour_discrete(name = "Immune\nprotein origin") +  
    ylab(expression(paste(italic("D. melanogaster"),  " protein length [aa]"))) + 
    xlab(expression(paste(italic("O. bicornis"), " protein length [aa]")))) 

(dm_stats2ob_stats_plot_log10 <- ggplot(dm_stats2ob_stats) + 
  geom_point(aes(x = Length_Ob, y = Length_Dm, colour = DiffNumDomains)) + 
  scale_colour_viridis_c(name = "Domain \ncount \ndifference", limits = c(-6, 6)) + 
  scale_x_log10() + 
  scale_y_log10() + ylab("D.melanogaster protein length (log10-transformed)") + 
  xlab("O.bicornis protein length (log10-transformed)"))

(bt_stats2ob_stats_plot <- ggplot(bt_stats2ob_stats) + 
  geom_point(aes(x = Length_Ob, y = Length_Bt, colour = DiffNumDomains)) + 
  scale_colour_viridis_c(name = "Domain \ncount \ndifference", limits = c(-6, 6)) +  
  ylab("B.terrestris protein length") + 
  xlab("O.bicornis protein length") + 
  xlim(0, max(c(bt_stats2ob_stats$Length_Bt, 
                bt_stats2ob_stats$Length_Ob) * 1.05)) + 
  ylim(0, max(c(bt_stats2ob_stats$Length_Bt, 
                bt_stats2ob_stats$Length_Ob) * 1.05)))

(bt_stats2ob_stats_plot_log10 <- ggplot(bt_stats2ob_stats) + 
  geom_point(aes(x = Length_Ob, y = Length_Bt, colour = DiffNumDomains)) + 
  scale_colour_viridis_c(name = "Domain \ncount \ndifference", limits = c(-6, 6)) + 
  scale_x_log10() + 
  scale_y_log10() + ylab("B.terrestris protein length (log10-transformed)") + 
  xlab("O.bicornis protein length (log10-transformed)"))

ggsave("results/PRJNA285788/immune_gene_identification/plots/prot_lengths_and_domains/dm_stats_vs_bt_stats.png", 
       dm_stats2bt_stats_plot)
ggsave("results/PRJNA285788/immune_gene_identification/plots/prot_lengths_and_domains/dm_stats_vs_bt_stats_log10.png", 
       dm_stats2bt_stats_plot_log10)
ggsave("results/PRJNA285788/immune_gene_identification/plots/prot_lengths_and_domains/dm_stats_vs_ob_stats_v1.png", 
       dm_stats2ob_stats_plot_v1)
ggsave("results/PRJNA285788/immune_gene_identification/plots/prot_lengths_and_domains/dm_stats_vs_ob_stats_v2.png", 
       dm_stats2ob_stats_plot_v2)
ggsave("results/PRJNA285788/immune_gene_identification/plots/prot_lengths_and_domains/dm_stats_vs_ob_stats_log10.png", 
       dm_stats2ob_stats_plot_log10)
ggsave("results/PRJNA285788/immune_gene_identification/plots/prot_lengths_and_domains/bt_stats_vs_ob_stats.png", 
       bt_stats2ob_stats_plot)
ggsave("results/PRJNA285788/immune_gene_identification/plots/prot_lengths_and_domains/bt_stats_vs_ob_stats_log10.png", 
       bt_stats2ob_stats_plot_log10)

print(sum(!dm_stats2bt_stats$AllEqualDomains & dm_stats2bt_stats$DiffNumDomains == 0))
print(sum(!dm_stats2ob_stats$AllEqualDomains & dm_stats2ob_stats$DiffNumDomains == 0))
print(sum(!bt_stats2ob_stats$AllEqualDomains & bt_stats2ob_stats$DiffNumDomains == 0))

save(overlap_imm_genes_plot_no_title, dm_stats2ob_stats_plot_v1, dm_stats2ob_stats_plot_v2, 
     unique_imm_gene_bt_seq_ident_distr_plot, unique_imm_gene_dm_seq_ident_distr_plot, 
     file = "results/PRJNA285788/immune_gene_identification/plots/figure1_plots.RData")

# Distributions -------------------------------------------------------------------------------

# ggsave("results/PRJNA285788/immune_gene_identification/plots/prot_lengths_and_domains/bt_ob_dom_num_distr.png", 
#        bt_ob_dom_num_distr_plot)
# ggsave("results/PRJNA285788/immune_gene_identification/plots/prot_lengths_and_domains/dm_ob_dom_num_distr.png", 
#        dm_ob_dom_num_distr_plot)
# ggsave("results/PRJNA285788/immune_gene_identification/plots/prot_lengths_and_domains/dm_bt_dom_num_distr.png", 
#        dm_bt_dom_num_distr_plot)
# 
# 

# transform_orthodb_df <- function(orthodb_df){
#   new_df <- tibble(Dm_ProteinID = "1", Bi_GeneID = "1")
#   for(i in 1:dim(orthodb_df)[1]){
#     for(dm_id in orthodb_df[["Dm_IDs"]][[i]]){
#       for(bi_id in orthodb_df[["Bi_IDs"]][[i]]){
#         new_df <- new_df %>% add_row(Dm_ProteinID = dm_id, Bi_GeneID = bi_id)
#       }
#     }
#   }
#   return(new_df[2:dim(new_df)[1],])
# }
# 
# dm_gene2dm_protein_fb <- as_tibble(biomaRt::getBM(c("ensembl_peptide_id", "ensembl_gene_id"), 
#                                                   mart = dm_mart)) %>% 
#   dplyr::filter(ensembl_peptide_id != "") %>% 
#   dplyr::rename("Dm_ProteinID" = "ensembl_peptide_id", "Dm_GeneID" = "ensembl_gene_id")
# 
# 
# bi_dm_orthodb_orthl_raw <- read.delim("data/reference/annotation/OrthoDB_Bi_Dm_orthologues.txt", 
#                                   sep = " ", header = FALSE) %>% as_tibble() %>% 
#   dplyr::select(V2, V4) %>% group_by(V2) %>% summarise(V4 = paste(V4, collapse = ", ")) %>% 
#   dplyr::filter(str_detect(V4, "BIMP") & str_detect(V4, "FBpp")) %>% 
#   mutate(Dm_IDs = str_extract_all(V4, "FBpp[0-9]+"), Bi_IDs = str_extract_all(V4, "BIMP[0-9]+")) %>% 
#   dplyr::select(-V4)
# 
# bi_dm_orthodb_orthl <- transform_orthodb_df(bi_dm_orthodb_orthl_raw)
# bi_dm_orthodb_orthl <- bi_dm_orthodb_orthl %>% left_join(dm_gene2dm_protein_fb) %>% 
#   dplyr::select(Dm_GeneID, Bi_GeneID) %>% distinct() %>% 
#   dplyr::filter(Dm_GeneID != "", !is.na(Dm_GeneID), Bi_GeneID != "", !is.na(Bi_GeneID))
# 
# bi_dm_ensembl_orthl <- biomaRt::getBM(c("ensembl_gene_id", "bimpatiens_eg_homolog_ensembl_gene"), 
#                                       mart = dm_mart) %>% as_tibble() %>% 
#   dplyr::rename("Dm_GeneID" = "ensembl_gene_id", 
#                 "Bi_GeneID" = "bimpatiens_eg_homolog_ensembl_gene") %>% 
#   dplyr::filter(!is.na(Dm_GeneID), Dm_GeneID != "", !is.na(Bi_GeneID), Bi_GeneID != "") %>% 
#   dplyr::distinct()
# 
# bi_dm_orthodb_orthl_4euler <- bi_dm_orthodb_orthl %>% 
#   mutate(forEuler = paste(Dm_GeneID, Bi_GeneID, sep = "to")) %>% dplyr::select(forEuler) %>% 
#   unlist() %>% as.vector()
# 
# bi_dm_ensembl_orthl_4euler <- bi_dm_ensembl_orthl %>% 
#   mutate(forEuler = paste(Dm_GeneID, Bi_GeneID, sep = "to")) %>% dplyr::select(forEuler) %>% 
#   unlist() %>% as.vector()
# 
# bi_dm_ensembl_orthodb_overlap_plot <- plot(eulerr::euler(list(
#   ensembl = bi_dm_ensembl_orthl_4euler, orthodb = bi_dm_orthodb_orthl_4euler)), quantities = TRUE)
# ggsave("results/PRJNA285788/immune_gene_identification/plots/prot_lengths_and_domains/bi_dm_ensembl_orthodb_overlap.png", 
#        bi_dm_ensembl_orthodb_overlap_plot)
# 
# transform_orthodb_df_am <- function(orthodb_df){
#   new_df <- tibble(Dm_ProteinID = "1", Am_ProteinID = "1")
#   for(i in 1:dim(orthodb_df)[1]){
#     for(dm_id in orthodb_df[["Dm_IDs"]][[i]]){
#       for(am_id in orthodb_df[["Am_IDs"]][[i]]){
#         new_df <- new_df %>% add_row(Dm_ProteinID = dm_id, Am_ProteinID = am_id)
#       }
#     }
#   }
#   return(new_df[2:dim(new_df)[1],])
# }
# 
# am_dm_orthodb_orthl_raw <- read.delim("data/reference/annotation/OrthoDB_Am_Dm_orthologues.txt", 
#                                       sep = " ", header = FALSE) %>% as_tibble() %>% 
#   dplyr::select(V2, V4) %>% group_by(V2) %>% summarise(V4 = paste(V4, collapse = ", ")) %>% 
#   dplyr::filter(str_detect(V4, "XP") & str_detect(V4, "FBpp")) %>% 
#   mutate(Am_IDs = str_extract_all(V4, "XP_[0-9]+"), Dm_IDs = str_extract_all(V4, "FBpp[0-9]+")) %>% 
#   dplyr::select(-V4)
# 
# am_dm_orthodb_orthl <- transform_orthodb_df_am(am_dm_orthodb_orthl_raw)
# am_dm_orthodb_orthl_gene_lvl <- am_dm_orthodb_orthl %>% left_join(dm_gene2dm_protein_fb) %>% 
#   left_join(am_gene2am_protein) %>% 
#   dplyr::select(Dm_GeneID, Am_GeneID) %>% distinct() %>% 
#   dplyr::filter(Dm_GeneID != "", !is.na(Dm_GeneID), Am_GeneID != "", !is.na(Am_GeneID))

