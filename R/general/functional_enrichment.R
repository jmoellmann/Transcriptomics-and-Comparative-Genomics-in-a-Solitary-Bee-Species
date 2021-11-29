library(tidyverse)
library(gggenes)
library(biomaRt)
library(topGO)

project_dir <- "/run/media/T5/masters_project/01_main_analysis/"
hyphy_res_dir <- paste0(project_dir, "01_data/02_phylogenetic_analysis/hyphy_results/")

setwd(project_dir)
source("02_scripts/PRJNA285788/general/get_gsea_results.R")

genes <- readRDS("04_results/PRJNA285788/genes.RDS")
single_copy_ogs <- read_table("01_data/03_reference/Orthogroups_SingleCopyOrthologues.txt")[[1]]

# Functional enrichment -----------------------------------------------------------------------

ob_selected_genes <- genes$Osmia_bicornis %>% filter(episDivSelection) %>% arrange(desc(dnDs)) %>% 
   pull(geneID)
ob_sco_genes <- genes$Osmia_bicornis %>% filter(!is.na(dnDs)) %>% pull(geneID)
ob_genes_pvals <- genes$Osmia_bicornis %>% filter(!is.na(dnDs)) %>% pull(selectionPval)
attr(ob_genes_pvals, "names") <- ob_sco_genes

ob_genes_ranked <- ob_genes_pvals %>% sort() %>% names()

dm_gene2ob_gene <- read_csv("04_results/PRJNA285788/immune_gene_identification/ob_gene2dm_gene.csv") %>%
   dplyr::rename("gene_id" = "Dm_GeneID")

dm_mart <- biomaRt::useEnsemblGenomes(biomart = "metazoa_mart",
                                      dataset = "dmelanogaster_eg_gene")
ob_gene2dm_GO <- biomaRt::getBM(c("ensembl_gene_id", "go_id", "go_linkage_type"),
                                mart = dm_mart) %>%
   as_tibble() %>%
   dplyr::rename("gene_id" = "ensembl_gene_id", "GOID" = "go_id",
                 "evidence" = "go_linkage_type") %>%
   left_join(dm_gene2ob_gene) %>%
   dplyr::select(Ob_GeneID, GOID) %>%
   dplyr::rename("Dm_GOID" = "GOID") %>% 
   drop_na()

ob_gene2dm_GO_list <- ob_gene2dm_GO %>% dplyr::filter(Ob_GeneID != "", Dm_GOID != "") %>%
   distinct() %>% as.data.frame() %>% GO_mappings_to_list()

ob_topGO_selected_ks <- get_gsea_results(ob_genes_pvals, ob_gene2dm_GO_list,
                                         outputTopGOobject = FALSE, nodeSize = 10, alpha = 0)
ob_topGO_selected_ks_depleted <- ob_topGO_selected_ks %>% filter(Significant < Expected) %>% 
   filter(classicKS < 0.05 & weight01KS < 0.05) %>% arrange(weight01KS)
ob_topGO_selected_ks_enriched <- ob_topGO_selected_ks %>% filter(Significant > Expected) %>% 
   filter(classicKS < 0.05 & weight01KS < 0.05) %>% arrange(weight01KS)
