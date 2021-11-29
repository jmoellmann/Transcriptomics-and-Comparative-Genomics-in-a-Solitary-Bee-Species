library(rjson)
library(tidyverse)
library(DESeq2)


project_dir <- "/run/media/T5/masters_project/01_main_analysis/"
hyphy_res_dir <- paste0(project_dir, "01_data/02_phylogenetic_analysis/hyphy_results/")

setwd(project_dir)

get_genes <- function(species_name){
   abbrev <- paste0(str_sub(str_to_lower(str_split(species_name, "_")[[1]][1]), 1, 1), 
                   str_sub(str_split(species_name, "_")[[1]][2], 1, 1))
   low_div_genes <- read.table(
      paste0("01_data/03_reference/low_diversity_genes/", abbrev, 
             "_low_diversity_genes_orthofinder.txt"))[[1]]
   gene2prot <- read_tsv(paste0("01_data/03_reference/annotation/", abbrev, "_gene2protein.tsv"), 
                            col_names = c("geneID", "protID"))
   gene2og <- read_tsv(
      "04_results/PRJNA285788/immune_gene_identification/orthofinder_bees/N0_new.tsv") %>% 
      dplyr::select(OG, species_name) %>% 
      dplyr::rename(protID = species_name) %>% 
      separate_rows(protID, sep = ",") %>% 
      inner_join(gene2prot) %>% dplyr::select(OG, geneID) %>% distinct()
   gene2chr_start_end <- read_tsv(paste0("01_data/03_reference/annotation/", abbrev, "_gene2chr_start_end.tsv"),
                                     col_names = c("geneID", "chr", "strand", "start", "end"))
   gene_dn_ds_tbl <- proteins_dn_ds_tbl[[species_name]] %>% inner_join(gene2prot) %>% 
      dplyr::select(geneID, dnDs, selectionPval) %>% distinct()
   selected_genes <- gene2prot %>% filter(protID %in% selected_prots[[species_name]]) %>% 
      pull(geneID)
   genes <- gene2chr_start_end %>% left_join(gene2og) %>% left_join(gene_dn_ds_tbl) %>% 
      mutate(episDivSelection = ifelse(geneID %in% selected_genes, TRUE, FALSE)) %>% 
      mutate(singleCopyOrthl = ifelse(OG %in% single_copy_ogs, TRUE, FALSE)) %>% 
      mutate(lowDivGene = ifelse(geneID %in% low_div_genes, TRUE, FALSE)) %>% 
      arrange(chr, start)
   return(genes)
}

# Get selected proteins -----------------------------------------------------------------------

proteins <- list()
proteins_pvals <- list()
selected_prots <- list()
proteins_dn_ds <- list()
selection_pvals <- list()
species_list <- c("Apis_mellifera", "Bombus_terrestris", "Dufourea_novaeangliae", "Megalopta_genalis", 
             "Eufriesea_mexicana", "Ceratina_calcarata", "Habropoda_laboriosa", "Megachile_rotundata",
             "Nomia_melanderi", "Osmia_bicornis", "Osmia_lignaria")

for (json in list.files(hyphy_res_dir)) {
   results <- fromJSON(file = paste0(hyphy_res_dir, json))[["branch attributes"]][["0"]]
   for (species in species_list){
      species_id <- names(results)[grep(species, names(results))]
      tmp_prot_id <- results[[species_id]][["original name"]]
      prot_id <- str_remove_all(tmp_prot_id, paste0(species, "_"))
      proteins[[species]] <- c(proteins[[species]], prot_id)
      avg_dn_ds <- 0
      dn_ds <- results[[species_id]][["Rate Distributions"]]
      for (site_specific_omega in dn_ds){
         avg_dn_ds = avg_dn_ds + site_specific_omega[1] * site_specific_omega[2]
      }
      selection_pval <- results[[species_id]][["Uncorrected P-value"]]
      selection_pvals[[species]] <- c(selection_pvals[[species]], selection_pval)
      proteins_dn_ds[[species]] <- c(proteins_dn_ds[[species]], avg_dn_ds)
      proteins_pvals[[species]] <- c(proteins_pvals[[species]], results[[species_id]][["Corrected P-Value"]])
      if(results[[species_id]][["Corrected P-value"]] < 0.05){
         selected_prots[[species]] <- c(selected_prots[[species]], prot_id)
      }
   }
}

proteins_dn_ds_tbl <- mapply(function(x,y,z) tibble(protID = x, dnDs = y, selectionPval = z), 
                             proteins, proteins_dn_ds, selection_pvals, SIMPLIFY = FALSE)

# Low diversity genes -------------------------------------------------------------------------

bt_low_div_genes <- read.table("01_data/03_reference/low_diversity_genes/bt_low_diversity_genes.txt")[[1]]
am_low_div_genes <- read.table("01_data/03_reference/low_diversity_genes/am_low_diversity_genes.txt")[[1]]
ob_low_div_genes_orthofinder <- read.table(
   "01_data/03_reference/low_diversity_genes/ob_low_diversity_genes_orthofinder.txt")[[1]]

single_copy_ogs <- read_table("01_data/03_reference/Orthogroups_SingleCopyOrthologues.txt")[[1]]

# Osmia bicornis ------------------------------------------------------------------------------

ob_gene2prot <- read_tsv("01_data/03_reference/annotation/ob_gene2protein.tsv", 
                         col_names = c("geneID", "protID"))

ob_imm_genes <- read_csv("04_results/PRJNA285788/immune_gene_identification/immune_genes_annotated.csv")[["geneID"]]

ob_spec_genes <- read_tsv("04_results/PRJNA285788/immune_gene_identification/osmia_specific_genes.tsv")[["GeneID"]]

ob_gene_dn_ds_tbl <- proteins_dn_ds_tbl[["Osmia_bicornis"]] %>% inner_join(ob_gene2prot) %>% 
   dplyr::select(geneID, dnDs, selectionPval) %>% distinct()

ob_selected_genes <- ob_gene2prot %>% filter(protID %in% selected_prots[["Osmia_bicornis"]]) %>% 
   pull(geneID)

ob_gene2og <- read_tsv(
   "04_results/PRJNA285788/immune_gene_identification/orthofinder_bees/N0_new.tsv") %>% 
   dplyr::select(OG, Osmia_bicornis) %>% 
   dplyr::rename(protID = Osmia_bicornis) %>% 
   separate_rows(protID, sep = ",") %>% 
   mutate(protID = str_remove(protID, " ")) %>%
   inner_join(ob_gene2prot) %>% dplyr::select(OG, geneID) %>% distinct()

ob_gene2chr_start_end <- read_tsv("01_data/03_reference/annotation/ob_gene2chr_start_end.tsv",
                                  col_names = c("geneID", "chr", "strand", "start", "end"))


avg_gene_counts <- readRDS("04_results/PRJNA285788/sex_differences/diff_expr/STAR/dds_lrt.RDS") %>% 
   DESeq2::counts(normalized = TRUE) %>% as_tibble(rownames = "geneID") %>% 
   mutate(meanReadCount = rowMeans(dplyr::select(., starts_with("SRR")))) %>% dplyr::select(geneID, meanReadCount)

avg_gene_counts_pesticides <- readRDS("04_results/PRJNA285788/pesticide_exposure/diff_expr/STAR/dds_wald.RDS") %>% 
   DESeq2::counts(normalized = TRUE) %>% as_tibble(rownames = "geneID") %>% 
   mutate(meanReadCount = rowMeans(dplyr::select(., starts_with("SRR")))) %>% dplyr::select(geneID, meanReadCount)

ob_genes <- read_csv("04_results/PRJNA285788/sex_differences/diff_expr/STAR/DESeq2_results.csv") %>% 
   dplyr::select(geneID, geneDescription, log2FoldChange, padj, pvalue) %>% left_join(ob_gene2prot) %>% 
   left_join(ob_gene2og) %>% left_join(ob_gene_dn_ds_tbl) %>% 
   mutate(episDivSelection = ifelse(geneID %in% ob_selected_genes, TRUE, FALSE)) %>% 
   group_by(geneID) %>% summarise(geneID, geneDescription, OG, log2FoldChange, padj, pvalue, 
                                  episDivSelection = any(episDivSelection), dnDs, selectionPval) %>% 
   distinct() %>% mutate(singleCopyOrthl = ifelse(OG %in% single_copy_ogs, TRUE, FALSE)) %>%
   mutate(lowDivGene = ifelse(geneID %in% ob_low_div_genes_orthofinder, TRUE, FALSE)) %>%
   mutate(immGene = ifelse(geneID %in% ob_imm_genes, TRUE, FALSE)) %>% 
   mutate(osmiaSpecific = ifelse(geneID %in% ob_spec_genes, TRUE, FALSE)) %>% 
   left_join(ob_gene2chr_start_end) %>% left_join(avg_gene_counts) %>% arrange(padj)

ob_genes_thiacloprid <- read_csv("04_results/PRJNA285788/pesticide_exposure/diff_expr/STAR/DESeq2_results_thiacloprid.csv") %>% 
   dplyr::select(geneID, geneDescription, log2FoldChange, padj, pvalue) %>% left_join(ob_gene2prot) %>% 
   left_join(ob_gene2og) %>% left_join(ob_gene_dn_ds_tbl) %>% 
   mutate(episDivSelection = ifelse(geneID %in% ob_selected_genes, TRUE, FALSE)) %>% 
   group_by(geneID) %>% summarise(geneID, geneDescription, OG, log2FoldChange, padj, pvalue, 
                                  episDivSelection = any(episDivSelection), dnDs, selectionPval) %>% 
   distinct() %>% mutate(singleCopyOrthl = ifelse(OG %in% single_copy_ogs, TRUE, FALSE)) %>%
   mutate(lowDivGene = ifelse(geneID %in% ob_low_div_genes_orthofinder, TRUE, FALSE)) %>%
   mutate(immGene = ifelse(geneID %in% ob_imm_genes, TRUE, FALSE)) %>% 
   mutate(osmiaSpecific = ifelse(geneID %in% ob_spec_genes, TRUE, FALSE)) %>% 
   left_join(ob_gene2chr_start_end) %>% left_join(avg_gene_counts_pesticides) %>% arrange(padj)

ob_genes_imidacloprid <- read_csv("04_results/PRJNA285788/pesticide_exposure/diff_expr/STAR/DESeq2_results_imidacloprid.csv") %>% 
   dplyr::select(geneID, geneDescription, log2FoldChange, padj, pvalue) %>% left_join(ob_gene2prot) %>% 
   left_join(ob_gene2og) %>% left_join(ob_gene_dn_ds_tbl) %>% 
   mutate(episDivSelection = ifelse(geneID %in% ob_selected_genes, TRUE, FALSE)) %>% 
   group_by(geneID) %>% summarise(geneID, geneDescription, OG, log2FoldChange, padj, pvalue, 
                                  episDivSelection = any(episDivSelection), dnDs, selectionPval) %>% 
   distinct() %>% mutate(singleCopyOrthl = ifelse(OG %in% single_copy_ogs, TRUE, FALSE)) %>%
   mutate(lowDivGene = ifelse(geneID %in% ob_low_div_genes_orthofinder, TRUE, FALSE)) %>%
   mutate(immGene = ifelse(geneID %in% ob_imm_genes, TRUE, FALSE)) %>% 
   mutate(osmiaSpecific = ifelse(geneID %in% ob_spec_genes, TRUE, FALSE)) %>% 
   left_join(ob_gene2chr_start_end) %>% left_join(avg_gene_counts_pesticides) %>% arrange(padj)

# Hypothesis testing for Osmia bicornis -------------------------------------------------------

male_biased_genes <- ob_genes %>% filter(singleCopyOrthl, padj < 0.05, log2FoldChange > 0)
female_biased_genes <- ob_genes %>% filter(singleCopyOrthl, padj < 0.05, log2FoldChange < 0)

prob_male_biased <- nrow(male_biased_genes) / (nrow(female_biased_genes) + nrow(male_biased_genes))

male_biased_selected_genes <- male_biased_genes %>% filter(episDivSelection)
female_biased_selected_genes <- female_biased_genes %>% filter(episDivSelection)

# Do we see more male biased genes under positive selection than female biased genes?
sex_biased_selection_test <- binom.test(x = nrow(male_biased_selected_genes), 
           n = nrow(male_biased_selected_genes) + nrow(female_biased_selected_genes),
           p = prob_male_biased)

# Do we see an enrichment of genes under positive selection among DEGs?
deg_enrichment_test_matrix <- matrix(data = c(nrow(ob_genes %>% filter(singleCopyOrthl)), 
                               nrow(ob_genes %>% filter(episDivSelection)),
                               nrow(ob_genes %>% filter(singleCopyOrthl, padj < 0.05)),
                               nrow(ob_genes %>% filter(episDivSelection, padj < 0.05))), nrow = 2)
deg_enrichment_test <- fisher.test(deg_enrichment_test_matrix)

# Do we see an enrichment of genes under positive selection among immune gene?
immune_enrichment_matrix <- matrix(data = c(nrow(ob_genes %>% filter(singleCopyOrthl)),
                                            nrow(ob_genes %>% filter(singleCopyOrthl & immGene)),
                                            nrow(ob_genes %>% filter(episDivSelection)),
                                            nrow(ob_genes %>% filter(episDivSelection & immGene))),
                                   nrow = 2)

immune_enrichment_test <- fisher.test(immune_enrichment_matrix)


# Apis mellifera ------------------------------------------------------------------------------

am_gene2prot <- read_tsv("01_data/03_reference/annotation/am_gene2protein_v2.tsv", 
                         col_names = c("geneID", "protID")) %>% drop_na()

am_prot2description <- read_tsv("01_data/03_reference/annotation/am_prot2description.tsv",
                                col_names = c("geneDescription", "protID"))

am_gene2description <- am_prot2description %>% inner_join(am_gene2prot) %>% 
   dplyr::select(geneID, geneDescription) %>% distinct()

am_gene_dn_ds_tbl <- proteins_dn_ds_tbl[["Apis_mellifera"]] %>% inner_join(am_gene2prot) %>% 
   dplyr::select(geneID, dnDs, selectionPval) %>% distinct()

am_selected_genes <- am_gene2prot %>% filter(protID %in% selected_prots[["Apis_mellifera"]]) %>% 
   pull(geneID)

am_gene2og <- read_tsv(
   "04_results/PRJNA285788/immune_gene_identification/orthofinder_bees/N0_new.tsv") %>% 
   dplyr::select(OG, Apis_mellifera) %>% 
   dplyr::rename(protID = Apis_mellifera) %>% 
   separate_rows(protID, sep = ",") %>% 
   inner_join(am_gene2prot) %>% dplyr::select(OG, geneID) %>% distinct()

am_gene2chr_start_end <- read_tsv("01_data/03_reference/annotation/am_gene2chr_start_end.tsv",
                                  col_names = c("geneID", "chr", "strand", "start", "end"))

am_genes <- results(readRDS("04_results/PRJNA338450/after_trimming/diff_expr/STAR/dds_lrt.RDS"))
#am_genes <- results(readRDS("04_results/PRJNA533756/after_trimming/diff_expr/STAR/dds_merged_sex.RDS"))
am_genes <- am_genes %>% as_tibble() %>% mutate(geneID = row.names(am_genes)) %>% 
   left_join(am_gene2prot) %>% left_join(am_gene2og) %>% left_join(am_gene2description) %>% 
   left_join(am_gene_dn_ds_tbl) %>% 
   mutate(episDivSelection = ifelse(geneID %in% am_selected_genes, TRUE, FALSE)) %>% 
   group_by(geneID) %>% summarise(geneID, geneDescription = geneDescription[1], OG, log2FoldChange, padj, 
                                  episDivSelection = episDivSelection, dnDs, selectionPval) %>% distinct() %>% 
   mutate(singleCopyOrthl = ifelse(OG %in% single_copy_ogs, TRUE, FALSE)) %>% 
   mutate(lowDivGene = ifelse(geneID %in% am_low_div_genes, TRUE, FALSE)) %>% 
   left_join(am_gene2chr_start_end) %>% arrange(padj)


# PCA on merged A.mellifera samples -----------------------------------------------------------

am_vst <- DESeq2::vst(readRDS("04_results/PRJNA533756/after_trimming/diff_expr/STAR/dds_merged_sex.RDS"))
DESeq2::plotPCA(am_vst, intgroup = "Group")

# Bombus terrestris ---------------------------------------------------------------------------

bt_gene2prot <- read_tsv("01_data/03_reference/annotation/bt_gene2protein.tsv", 
                         col_names = c("geneID", "protID"))

bt_gene2og <- read_tsv(
   "04_results/PRJNA285788/immune_gene_identification/orthofinder_bees/N0_new.tsv") %>% 
   dplyr::select(OG, Bombus_terrestris) %>% 
   dplyr::rename(protID = Bombus_terrestris) %>% 
   separate_rows(protID, sep = ",") %>% 
   inner_join(bt_gene2prot) %>% dplyr::select(OG, geneID) %>% distinct()

bt_gene2description <- read_tsv("01_data/03_reference/annotation/bt_gene2gene_description.tsv",
         col_names = c("geneID", "geneDescription"))


bt_gene_dn_ds_tbl <- proteins_dn_ds_tbl[["Bombus_terrestris"]] %>% inner_join(bt_gene2prot) %>% 
   dplyr::select(geneID, dnDs, selectionPval) %>% distinct()

bt_selected_genes <- bt_gene2prot %>% filter(protID %in% selected_prots[["Bombus_terrestris"]]) %>% 
   pull(geneID)

bt_gene2chr_start_end <- read_tsv("01_data/03_reference/annotation/bt_gene2chr_start_end.tsv",
                                  col_names = c("geneID", "chr", "strand", "start", "end"))

bt_genes <- results(readRDS("04_results/PRJEB9366/diff_expr/STAR/dds_lrt.RDS"))
bt_genes <- bt_genes %>% as_tibble() %>% mutate(geneID = row.names(bt_genes)) %>% 
   left_join(bt_gene2description) %>% left_join(bt_gene2og) %>% 
   left_join(bt_gene2prot) %>% left_join(bt_gene_dn_ds_tbl) %>% 
   mutate(episDivSelection = ifelse(geneID %in% bt_selected_genes, TRUE, FALSE)) %>% 
   group_by(geneID) %>% summarise(geneID, geneDescription = geneDescription[1], OG, log2FoldChange, padj, 
                                  episDivSelection = any(episDivSelection), dnDs, selectionPval) %>% distinct() %>% 
   mutate(singleCopyOrthl = ifelse(OG %in% single_copy_ogs, TRUE, FALSE)) %>% 
   mutate(lowDivGene = ifelse(geneID %in% bt_low_div_genes, TRUE, FALSE)) %>% 
   left_join(bt_gene2chr_start_end) %>% 
   arrange(padj)

bt_vst <- DESeq2::vst(readRDS("04_results/PRJEB9366/diff_expr/STAR/dds_wald.RDS"))
DESeq2::plotPCA(bt_vst, intgroup = "group")



# Other species -------------------------------------------------------------------------------

genes <- list()
species_names <- c("Megalopta_genalis", "Dufourea_novaeangliae", "Ceratina_calcarata", 
                   "Eufriesea_mexicana", "Habropoda_laboriosa", "Megachile_rotundata", 
                   "Nomia_melanderi", "Osmia_lignaria")
for(species_name in species_names){
   genes[[species_name]] <- get_genes(species_name)
}

genes[["Osmia_bicornis"]] <- ob_genes
genes[["Apis_mellifera"]] <- am_genes
genes[["Bombus_terrestris"]] <- bt_genes


# save ----------------------------------------------------------------------------------------

saveRDS(genes, "04_results/PRJNA285788/genes.RDS")

saveRDS(ob_genes_thiacloprid, "04_results/PRJNA285788/thiacloprid_genes.RDS")
saveRDS(ob_genes_imidacloprid, "04_results/PRJNA285788/imidacloprid_genes.RDS")

