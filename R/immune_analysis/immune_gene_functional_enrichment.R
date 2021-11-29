library(tidyverse)
library(topGO)

project_dir <- "/media/jannik/Samsung_T5/masters_project/ApisTranscriptomics/"
setwd(project_dir)

source("scripts/PRJNA285788/general/get_gsea_results.R")

res_dir <- "results/PRJNA285788/immune_gene_identification/funct_enrichment"
dir.create(res_dir, 
           showWarnings = FALSE, recursive = TRUE)

# Load immune genes ---------------------------------------------------------------------------

bt_imm_genes <- read.table("results/PRJNA285788/immune_gene_identification/immune_genes/bt_imm_genes.txt", 
           header = FALSE) %>% unlist() %>% as.vector()
dm_imm_genes <- read.table("results/PRJNA285788/immune_gene_identification/immune_genes/dm_imm_genes.txt", 
           header = FALSE) %>% unlist() %>% as.vector()
am_imm_genes <- read.table("results/PRJNA285788/immune_gene_identification/immune_genes/am_imm_genes.txt", 
           header = FALSE) %>% unlist() %>% as.vector()
dm_am_imm_genes <- read.table("results/PRJNA285788/immune_gene_identification/immune_genes/dm_am_imm_genes.txt", 
           header = FALSE) %>% unlist() %>% as.vector()
ob_bt_imm_genes <- read.table("results/PRJNA285788/immune_gene_identification/immune_genes/ob_bt_imm_genes.txt", 
           header = FALSE) %>% unlist() %>% as.vector()
ob_dm_imm_genes <- read.table("results/PRJNA285788/immune_gene_identification/immune_genes/ob_dm_imm_genes.txt", 
           header = FALSE) %>% unlist() %>% as.vector()
ob_am_imm_genes <- read.table("results/PRJNA285788/immune_gene_identification/immune_genes/ob_am_imm_genes.txt", 
           header = FALSE) %>% unlist() %>% as.vector()
bt_dm_imm_genes <- read.table("results/PRJNA285788/immune_gene_identification/immune_genes/bt_dm_imm_genes.txt", 
           header = FALSE) %>% unlist() %>% as.vector()
dm_bt_imm_genes <- read.table("results/PRJNA285788/immune_gene_identification/immune_genes/dm_bt_imm_genes.txt", 
           header = FALSE) %>% unlist() %>% as.vector()
ob_merged_imm_genes <- read.table("results/PRJNA285788/immune_gene_identification/immune_genes/ob_merged_imm_genes.txt", 
           header = FALSE) %>% unlist() %>% as.vector()


# Load Ob Gene to Dm/Bt GO terms tables -------------------------------------------------------

ob_gene2dm_GO <- read_csv("results/PRJNA285788/immune_gene_identification/funct_enrichment/ob_gene2dm_GO.csv")
ob_gene2bt_GO <- read_csv("results/PRJNA285788/immune_gene_identification/funct_enrichment/ob_gene2bt_GO.csv")

ob_gene2dm_GO_nona <- ob_gene2dm_GO %>% dplyr::filter(Ob_GeneID != "", Dm_GOID != "") %>% 
  distinct()
ob_gene2dm_GO_list <- GO_mappings_to_list(as.data.frame(ob_gene2dm_GO_nona))

ob_gene2bt_GO_nona <- ob_gene2bt_GO %>% dplyr::filter(Ob_GeneID != "", Bt_GOID != "") %>% 
  distinct()
ob_gene2bt_GO_list <- GO_mappings_to_list(as.data.frame(ob_gene2bt_GO_nona))

# bt_genes <- unique(bt_gene2bt_GO$Bt_GeneID)
# bt_gene_universe <- factor(as.integer(bt_genes %in% bt_imm_genes)) 
# names(bt_gene_universe) <- bt_genes
# bt_dm_gene_universe <- factor(as.integer(bt_genes %in% bt_dm_imm_genes)) 
# names(bt_dm_gene_universe) <- bt_genes
# 
# dm_genes <- unique(dm_gene2dm_GO$Dm_GeneID)
# dm_gene_universe <- factor(as.integer(dm_genes %in% dm_imm_genes)) 
# names(dm_gene_universe) <- dm_genes
# dm_bt_gene_universe <- factor(as.integer(dm_genes %in% dm_bt_imm_genes)) 
# names(dm_bt_gene_universe) <- dm_genes
# dm_am_gene_universe <- factor(as.integer(dm_genes %in% dm_am_imm_genes)) 
# names(dm_am_gene_universe) <- dm_genes

ob_genes <- unique(ob_gene2dm_GO$Ob_GeneID)
ob_dm_gene_universe <- factor(as.integer(ob_genes %in% ob_dm_imm_genes)) 
names(ob_dm_gene_universe) <- ob_genes
ob_bt_gene_universe <- factor(as.integer(ob_genes %in% ob_bt_imm_genes))
names(ob_bt_gene_universe) <- ob_genes
ob_am_gene_universe <- factor(as.integer(ob_genes %in% ob_am_imm_genes))
names(ob_am_gene_universe) <- ob_genes
ob_merged_gene_universe <- factor(as.integer(ob_genes %in% ob_merged_imm_genes))
names(ob_merged_gene_universe) <- ob_genes
ob_dm_bt_overlap_gene_universe <- factor(as.integer(ob_genes %in% ob_bt_imm_genes & 
                                                      ob_genes %in% ob_dm_imm_genes))
names(ob_dm_bt_overlap_gene_universe) <- ob_genes
ob_dm_unique_gene_universe <- factor(as.integer(ob_genes %in% ob_dm_imm_genes &
                                                  !(ob_genes %in% ob_bt_imm_genes)))
names(ob_dm_unique_gene_universe) <- ob_genes
ob_bt_unique_gene_universe <- factor(as.integer(ob_genes %in% ob_bt_imm_genes &
                                                  !(ob_genes %in% ob_dm_imm_genes)))
names(ob_bt_unique_gene_universe) <- ob_genes

ob_am_unique_gene_universe <- factor(as.integer(ob_genes %in% ob_am_imm_genes &
                                                  !(ob_genes %in% ob_dm_imm_genes)))
names(ob_am_unique_gene_universe) <- ob_genes

for(ONTOLOGY in c("BP", "MF", "CC")){
  
  dir.create(paste0(res_dir, "/", ONTOLOGY), showWarnings = FALSE, recursive = TRUE)
  
  # bt_imm_genes_topGO <- get_gsea_results(stats = bt_gene_universe, 
  #                                        gene2GO = bt_gene2bt_GO_list, mode = "fisher")
  # bt_dm_imm_genes_topGO <- get_gsea_results(stats = bt_dm_gene_universe, 
  #                                           gene2GO = bt_gene2bt_GO_list, mode = "fisher")
  # write_csv(bt_imm_genes_topGO, paste0(res_dir, "/", ONTOLOGY, "/bt_imm_genes_topGO.csv"))
  # write_csv(bt_dm_imm_genes_topGO, paste0(res_dir, "/", ONTOLOGY, "/bt_dm_imm_genes_topGO.csv"))
  # 
  # dm_imm_genes_topGO <- get_gsea_results(stats = dm_gene_universe, 
  #                                        gene2GO = dm_gene2dm_GO_list, mode = "fisher")
  # dm_bt_imm_genes_topGO <- get_gsea_results(stats = dm_bt_gene_universe, 
  #                                           gene2GO = dm_gene2dm_GO_list, mode = "fisher")
  # dm_am_imm_genes_topGO <- get_gsea_results(stats = dm_am_gene_universe, 
  #                                           gene2GO = dm_gene2dm_GO_list, mode = "fisher")
  # write_csv(dm_imm_genes_topGO, 
  #           paste0(res_dir, "/", ONTOLOGY, "/dm_imm_genes_topGO.csv"))
  # write_csv(dm_bt_imm_genes_topGO, 
  #           paste0(res_dir, "/", ONTOLOGY, "/dm_bt_imm_genes_topGO.csv"))
  # write_csv(dm_am_imm_genes_topGO, 
  #           paste0(res_dir, "/", ONTOLOGY, "/dm_am_imm_genes_topGO.csv"))
  
  ob_dm_imm_genes_topGO <- get_gsea_results(stats = ob_dm_gene_universe, 
                                            gene2GO = ob_gene2dm_GO_list, mode = "fisher",
                                            ontology = ONTOLOGY)
  ob_bt_imm_genes_topGO <- get_gsea_results(stats = ob_bt_gene_universe, 
                                            gene2GO = ob_gene2dm_GO_list, mode = "fisher",
                                            ontology = ONTOLOGY)
  ob_am_imm_genes_topGO <- get_gsea_results(stats = ob_am_gene_universe, 
                                            gene2GO = ob_gene2dm_GO_list, mode = "fisher",
                                            ontology = ONTOLOGY)
  ob_merged_imm_genes_topGO <- get_gsea_results(stats = ob_merged_gene_universe, 
                                                gene2GO = ob_gene2dm_GO_list, mode = "fisher",
                                                ontology = ONTOLOGY)
  ob_dm_bt_overlap_imm_genes_topGO <- get_gsea_results(stats = ob_dm_bt_overlap_gene_universe, 
                                                       gene2GO = ob_gene2dm_GO_list, mode = "fisher",
                                                       ontology = ONTOLOGY)
  ob_bt_unique_imm_genes_topGO <- get_gsea_results(stats = ob_bt_unique_gene_universe, 
                                                   gene2GO = ob_gene2dm_GO_list, mode = "fisher",
                                                   ontology = ONTOLOGY)
  ob_dm_unique_imm_genes_topGO <- get_gsea_results(stats = ob_dm_unique_gene_universe, 
                                                   gene2GO = ob_gene2dm_GO_list, mode = "fisher",
                                                   ontology = ONTOLOGY)
  ob_am_unique_imm_genes_topGO <- get_gsea_results(stats = ob_am_unique_gene_universe, 
                                                   gene2GO = ob_gene2dm_GO_list, mode = "fisher",
                                                   ontology = ONTOLOGY)
  
  write_csv(ob_dm_imm_genes_topGO, 
            paste0(res_dir, "/", "/ob_dm_imm_genes_topGO.csv"))
  write_csv(ob_bt_imm_genes_topGO, 
            paste0(res_dir, "/", ONTOLOGY, "/ob_bt_imm_genes_topGO.csv"))
  write_csv(ob_am_imm_genes_topGO, 
            paste0(res_dir, "/", ONTOLOGY, "/ob_am_imm_genes_topGO.csv"))
  write_csv(ob_dm_bt_overlap_imm_genes_topGO, 
            paste0(res_dir, "/", ONTOLOGY, "/ob_dm_bt_overlap_imm_genes_topGO.csv"))
  write_csv(ob_bt_unique_imm_genes_topGO, 
            paste0(res_dir, "/", ONTOLOGY, "/ob_bt_unique_imm_genes_topGO.csv"))
  write_csv(ob_dm_unique_imm_genes_topGO, 
            paste0(res_dir, "/", ONTOLOGY, "/ob_dm_unique_imm_genes_topGO.csv"))
  write_csv(ob_am_unique_imm_genes_topGO, 
            paste0(res_dir, "/", ONTOLOGY, "/ob_am_unique_imm_genes_topGO.csv"))
  write_csv(ob_merged_imm_genes_topGO, 
            paste0(res_dir, "/", ONTOLOGY, "/ob_merged_imm_genes_topGO.csv"))
}
