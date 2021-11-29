library(DESeq2)
library(ggplot2)
library(topGO)
library(dplyr)
library(stringr)
library(readr)
library(biomaRt)

project_dir <- "/run/media/T5/masters_project/01_main_analysis/"
setwd(project_dir)

source("02_scripts/PRJNA285788/general/get_gsea_results.R")
# Read in immune gene orthologues -------------------------------------------------------------

immune_genes <- read.table(
  "04_results/PRJNA285788/immune_gene_identification/immune_genes/ob_merged_imm_genes.txt")[["V1"]]

# Load Ob Gene to Dm/Bt GO terms tables -------------------------------------------------------

ob_gene2dm_GO <- read_csv(
   "04_results/PRJNA285788/immune_gene_identification/funct_enrichment/ob_gene2dm_GO.csv")
# ob_gene2bt_GO <- read_csv(
#    "04_results/PRJNA285788/immune_gene_identification/funct_enrichment/ob_gene2bt_GO.csv")

ob_gene2dm_GO_nona <- ob_gene2dm_GO %>% dplyr::filter(Ob_GeneID != "", Dm_GOID != "") %>% 
  distinct()
ob_gene2dm_GO_list <- GO_mappings_to_list(as.data.frame(ob_gene2dm_GO_nona))

# ob_gene2bt_GO_nona <- ob_gene2bt_GO %>% dplyr::filter(Ob_GeneID != "", Bt_GOID != "") %>% 
#   distinct()
# ob_gene2bt_GO_list <- GO_mappings_to_list(as.data.frame(ob_gene2bt_GO_nona))

res <- DESeq2::results(readRDS(paste0(
  "04_results/PRJNA285788/sex_differences/diff_expr/STAR/dds_lrt.RDS")))

# Run GO enrichment analysis ------------------------------------------------------------------
for(ontology in c("BP", "MF", "CC")){
  
  dir.create(paste0("04_results/PRJNA285788/sex_differences/diff_expr/STAR/funct_enrich/", ontology), 
             showWarnings = FALSE, recursive = TRUE)
  dir.create(paste0("04_results/PRJNA285788/sex_differences/diff_expr/kallisto/funct_enrich/", ontology), 
             showWarnings = FALSE, recursive = TRUE)
  
  
  stats <- res$pvalue
  names(stats) <- row.names(res)
  stats <- stats[!is.na(stats)]
  
  imm_gene_stats <- stats[names(stats) %in% immune_genes]
  
  gsea_ns20 <- get_gsea_results(stats, ob_gene2dm_GO_list, nodeSize = 20, 
                                ontology = ontology, mode = "both", alpha = 0.05)
  # gsea_ns50 <- get_gsea_results(stats, ob_gene2dm_GO_list, nodeSize = 50, 
  #                               ontology = ontology)
  # gsea_ns100 <- get_gsea_results(stats, ob_gene2dm_GO_list, nodeSize = 100, 
  #                                ontology = ontology)
  # gsea_ns150 <- get_gsea_results(stats, ob_gene2dm_GO_list, nodeSize = 150, 
  #                                ontology = ontology)
  # 
  # gsea_bt_ns10 <- get_gsea_results(stats, ob_gene2bt_GO_list, nodeSize = 10, 
  #                                  ontology = ontology)
  # gsea_bt_ns20 <- get_gsea_results(stats, ob_gene2bt_GO_list, nodeSize = 20, 
  #                                  ontology = ontology)
  # 
  # gsea_imm_gene_ns10 <- get_gsea_results(imm_gene_stats, ob_gene2dm_GO_list, nodeSize = 10, 
  #                                        ontology = ontology)
  gsea_imm_gene_ns20 <- get_gsea_results(imm_gene_stats, ob_gene2dm_GO_list, nodeSize = 20, 
                                         ontology = ontology, mode = "both", alpha = 0.05)
  # gsea_imm_gene_ns50 <- get_gsea_results(imm_gene_stats, ob_gene2dm_GO_list, nodeSize = 50, 
  #                                        ontology = ontology)
  # gsea_imm_gene_ns100 <- get_gsea_results(imm_gene_stats, ob_gene2dm_GO_list, nodeSize = 100, 
  #                                         ontology = ontology)
  # gsea_imm_gene_ns150 <- get_gsea_results(imm_gene_stats, ob_gene2dm_GO_list, nodeSize = 150, 
  #                                         ontology = ontology)
  # 
  # gsea_imm_gene_bt_ns10 <- get_gsea_results(imm_gene_stats, ob_gene2bt_GO_list, nodeSize = 10, 
  #                                           ontology = ontology)
  

  write_tsv(gsea_ns20, paste0("04_results/PRJNA285788/sex_differences/diff_expr/STAR/funct_enrich/", 
                              ontology, "/gsea_dm_ns20.tsv"))
  # write_tsv(gsea_ns50, paste0("results/PRJNA285788/sex_differences/diff_expr/STAR/funct_enrich/", 
  #                             ontology, "/gsea_dm_ns50.tsv"))
  # write_tsv(gsea_ns100, paste0("results/PRJNA285788/sex_differences/diff_expr/STAR/funct_enrich/",
  #                              ontology, "/gsea_dm_ns100.tsv"))
  # write_tsv(gsea_ns150, paste0("results/PRJNA285788/sex_differences/diff_expr/STAR/funct_enrich/",
  #                              ontology, "/gsea_dm_ns150.tsv"))
  # 
  # write_tsv(gsea_bt_ns10, paste0("results/PRJNA285788/sex_differences/diff_expr/STAR/funct_enrich/",
  #                                ontology, "/gsea_bt_ns10.tsv"))
  # write_tsv(gsea_bt_ns20, paste0("results/PRJNA285788/sex_differences/diff_expr/STAR/funct_enrich/", 
  #                                ontology, "/gsea_bt_ns20.tsv"))
  
  # write_tsv(gsea_imm_gene_ns10, paste0("results/PRJNA285788/sex_differences/diff_expr/STAR/funct_enrich/", 
  #                                      ontology, "/gsea_imm_gene_ns10.tsv"))
  write_tsv(gsea_imm_gene_ns20, paste0("04_results/PRJNA285788/sex_differences/diff_expr/STAR/funct_enrich/", 
                                       ontology, "/gsea_imm_gene_ns20.tsv"))
  # write_tsv(gsea_imm_gene_ns50, paste0("results/PRJNA285788/sex_differences/diff_expr/STAR/funct_enrich/", 
  #                                      ontology, "/gsea_imm_gene_ns50.tsv"))
  # write_tsv(gsea_imm_gene_ns100, paste0("results/PRJNA285788/sex_differences/diff_expr/STAR/funct_enrich/", 
  #                                       ontology, "/gsea_imm_gene_ns100.tsv"))
  # write_tsv(gsea_imm_gene_ns150, paste0("results/PRJNA285788/sex_differences/diff_expr/STAR/funct_enrich/", 
  #                                       ontology, "/gsea_imm_gene_ns150.tsv"))
  # 
  # write_tsv(gsea_imm_gene_bt_ns10, paste0("results/PRJNA285788/sex_differences/diff_expr/STAR/funct_enrich/", 
  #                                         ontology, "/gsea_imm_gene_bt_ns10.tsv"))
  
  stats2 <- res$log2FoldChange
  names(stats2) <- row.names(res)
  stats2 <- stats2[!is.na(stats2)]
  
  imm_gene_stats2 <- stats2[names(stats2) %in% immune_genes]
  
  gsea_ns20 <- get_gsea_results(stats2, ob_gene2dm_GO_list, nodeSize = 20, 
                                ontology = ontology, mode = "ks", alpha = 0)
  # gsea_ns50 <- get_gsea_results(stats2, ob_gene2dm_GO_list, nodeSize = 50, 
  #                               ontology = ontology, alpha = 0)
  # gsea_ns100 <- get_gsea_results(stats2, ob_gene2dm_GO_list, nodeSize = 100, 
  #                                ontology = ontology, alpha = 0)
  # 
  # gsea_bt_ns10 <- get_gsea_results(stats2, ob_gene2bt_GO_list, nodeSize = 10, 
  #                                  ontology = ontology, alpha = 0)
  # gsea_bt_ns20 <- get_gsea_results(stats2, ob_gene2bt_GO_list, nodeSize = 20, 
  #                                  ontology = ontology, alpha = 0)
  
  # gsea_imm_gene_ns10 <- get_gsea_results(imm_gene_stats2, ob_gene2dm_GO_list, nodeSize = 10, 
  #                                        ontology = ontology, alpha = 0)
  gsea_imm_gene_ns20 <- get_gsea_results(imm_gene_stats2, ob_gene2dm_GO_list, nodeSize = 20, 
                                         ontology = ontology, mode = "ks", alpha = 0)
  # gsea_imm_gene_ns50 <- get_gsea_results(imm_gene_stats2, ob_gene2dm_GO_list, nodeSize = 50, 
  #                                        ontology = ontology, alpha = 0)
  # gsea_imm_gene_ns100 <- get_gsea_results(imm_gene_stats2, ob_gene2dm_GO_list, nodeSize = 100, 
  #                                         ontology = ontology, alpha = 0)
  # 
  # gsea_imm_gene_bt_ns10 <- get_gsea_results(imm_gene_stats2, ob_gene2bt_GO_list, nodeSize = 10, 
  #                                           ontology = ontology, alpha = 0)
  # 
  # 
  write_tsv(gsea_ns20, paste0("04_results/PRJNA285788/sex_differences/diff_expr/STAR/funct_enrich/", 
                              ontology, "/gsea_dm_ns20_lfc.tsv"))
  # write_tsv(gsea_ns50, paste0("results/PRJNA285788/sex_differences/diff_expr/STAR/funct_enrich/", 
  #                             ontology, "/gsea_dm_ns50_lfc.tsv"))
  # write_tsv(gsea_ns100, paste0("results/PRJNA285788/sex_differences/diff_expr/STAR/funct_enrich/", 
  #                              ontology, "/gsea_dm_ns100_lfc.tsv"))
  
  # write_tsv(gsea_bt_ns10, paste0("results/PRJNA285788/sex_differences/diff_expr/STAR/funct_enrich/", 
  #                                ontology, "/gsea_bt_ns10_lfc.tsv"))
  # write_tsv(gsea_bt_ns20, paste0("results/PRJNA285788/sex_differences/diff_expr/STAR/funct_enrich/", 
  #                                ontology, "/gsea_bt_ns20_lfc.tsv"))
  # 
  # write_tsv(gsea_imm_gene_ns10, paste0("results/PRJNA285788/sex_differences/diff_expr/STAR/funct_enrich/", 
  #                                      ontology, "/gsea_imm_gene_ns10_lfc.tsv"))
  write_tsv(gsea_imm_gene_ns20, paste0("04_results/PRJNA285788/sex_differences/diff_expr/STAR/funct_enrich/", 
                                       ontology, "/gsea_imm_gene_ns20_lfc.tsv"))
  # write_tsv(gsea_imm_gene_ns50, paste0("results/PRJNA285788/sex_differences/diff_expr/STAR/funct_enrich/", 
  #                                      ontology, "/gsea_imm_gene_ns50_lfc.tsv"))
  # write_tsv(gsea_imm_gene_ns100, paste0("results/PRJNA285788/sex_differences/diff_expr/STAR/funct_enrich/", 
  #                                       ontology, "/gsea_imm_gene_ns100_lfc.tsv"))
  # 
  # write_tsv(gsea_imm_gene_bt_ns10, paste0("results/PRJNA285788/sex_differences/diff_expr/STAR/funct_enrich/", 
  #                                         ontology, "/gsea_imm_gene_bt_ns10_lfc.tsv"))
  
  stats3 <- res$log2FoldChange * (-1)
  names(stats3) <- row.names(res)
  stats3 <- stats3[!is.na(stats3)]
  
  imm_gene_stats3 <- stats3[names(stats3) %in% immune_genes]
  
  gsea_ns20 <- get_gsea_results(stats3, ob_gene2dm_GO_list, nodeSize = 20, 
                                ontology = ontology, mode = "ks", alpha = 0)
  # gsea_ns50 <- get_gsea_results(stats3, ob_gene2dm_GO_list, nodeSize = 50, 
  #                               ontology = ontology, alpha = 0)
  # gsea_ns100 <- get_gsea_results(stats3, ob_gene2dm_GO_list, nodeSize = 100, 
  #                                ontology = ontology, alpha = 0)
  # # 
  # gsea_bt_ns10 <- get_gsea_results(stats3, ob_gene2bt_GO_list, nodeSize = 10, 
  #                                  ontology = ontology, alpha = 0)
  # gsea_bt_ns20 <- get_gsea_results(stats3, ob_gene2bt_GO_list, nodeSize = 20, 
  #                                  ontology = ontology, alpha = 0)
  # 
  # gsea_imm_gene_ns10 <- get_gsea_results(imm_gene_stats3, ob_gene2dm_GO_list, nodeSize = 10, 
  #                                        ontology = ontology, alpha = 0)
  gsea_imm_gene_ns20 <- get_gsea_results(imm_gene_stats3, ob_gene2dm_GO_list, nodeSize = 20, 
                                         ontology = ontology, mode = "ks", alpha = 0)
  # gsea_imm_gene_ns50 <- get_gsea_results(imm_gene_stats3, ob_gene2dm_GO_list, nodeSize = 50, 
  #                                        ontology = ontology, alpha = 0)
  # gsea_imm_gene_ns100 <- get_gsea_results(imm_gene_stats3, ob_gene2dm_GO_list, nodeSize = 100, 
  #                                         ontology = ontology, alpha = 0)
  # 
  # gsea_imm_gene_bt_ns10 <- get_gsea_results(imm_gene_stats3, ob_gene2bt_GO_list, nodeSize = 10, 
  #                                           ontology = ontology, alpha = 0)
  # 
  
  write_tsv(gsea_ns20, paste0("04_results/PRJNA285788/sex_differences/diff_expr/STAR/funct_enrich/", 
                              ontology, "/gsea_dm_ns20_lfc_neg.tsv"))
  # write_tsv(gsea_ns50, paste0("results/PRJNA285788/sex_differences/diff_expr/STAR/funct_enrich/", 
  #                             ontology, "/gsea_dm_ns50_lfc_neg.tsv"))
  # write_tsv(gsea_ns100, paste0("results/PRJNA285788/sex_differences/diff_expr/STAR/funct_enrich/", 
  #                              ontology, "/gsea_dm_ns100_lfc_neg.tsv"))
  # 
  # write_tsv(gsea_bt_ns10, paste0("results/PRJNA285788/sex_differences/diff_expr/STAR/funct_enrich/", 
  #                                ontology, "/gsea_bt_ns10_lfc_neg.tsv"))
  # write_tsv(gsea_bt_ns20, paste0("results/PRJNA285788/sex_differences/diff_expr/STAR/funct_enrich/", 
  #                                ontology, "/gsea_bt_ns20_lfc_neg.tsv"))
  # 
  # write_tsv(gsea_imm_gene_ns10, paste0("results/PRJNA285788/sex_differences/diff_expr/STAR/funct_enrich/", 
  #                                      ontology, "/gsea_imm_gene_ns10_lfc_neg.tsv"))
  write_tsv(gsea_imm_gene_ns20, paste0("04_results/PRJNA285788/sex_differences/diff_expr/STAR/funct_enrich/", 
                                       ontology, "/gsea_imm_gene_ns20_lfc_neg.tsv"))
  # write_tsv(gsea_imm_gene_ns50, paste0("results/PRJNA285788/sex_differences/diff_expr/STAR/funct_enrich/", 
  #                                      ontology, "/gsea_imm_gene_ns50_lfc_neg.tsv"))
  # write_tsv(gsea_imm_gene_ns100, paste0("results/PRJNA285788/sex_differences/diff_expr/STAR/funct_enrich/", 
  #                                       ontology, "/gsea_imm_gene_ns100_lfc_neg.tsv"))
  # 
  # write_tsv(gsea_imm_gene_bt_ns10, paste0("results/PRJNA285788/sex_differences/diff_expr/STAR/funct_enrich/", 
  #                                         ontology, "/gsea_imm_gene_bt_ns10_lfc_neg.tsv"))
}
