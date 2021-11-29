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
ob_selected_genes_pvals <- genes$Osmia_bicornis %>% filter(!is.na(dnDs)) %>% pull(selectionPval)
attr(ob_selected_genes_pvals, "names") <- ob_sco_genes

ob_genes_ranked <- ob_selected_genes_pvals %>% sort() %>% names()

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


# BP ------------------------------------------------------------------------------------------

ob_topGO_selected_ks_bp <- get_gsea_results(ob_selected_genes_pvals, ob_gene2dm_GO_list, mode = "both",
                                         outputTopGOobject = FALSE, nodeSize = 10, ontology = "BP")
# ob_topGO_selected_ks_depleted <- ob_topGO_selected_ks %>% filter(Significant < Expected) %>% 
#    filter(classicFisher < 0.05) %>% arrange(classicFisher)
ob_topGO_selected_ks_bp_enriched <- ob_topGO_selected_ks_bp %>% filter(Significant > Expected) %>% 
   filter(classicFisher < 0.05) %>% arrange(classicFisher) %>% mutate(ontology = "BP")


# MF ------------------------------------------------------------------------------------------

ob_topGO_selected_ks_mf <- get_gsea_results(ob_selected_genes_pvals, ob_gene2dm_GO_list, mode = "both",
                                            outputTopGOobject = FALSE, nodeSize = 10, ontology = "MF")
# ob_topGO_selected_ks_depleted <- ob_topGO_selected_ks %>% filter(Significant < Expected) %>% 
#    filter(classicFisher < 0.05) %>% arrange(classicFisher)
ob_topGO_selected_ks_mf_enriched <- ob_topGO_selected_ks_mf %>% filter(Significant > Expected) %>% 
   filter(classicFisher < 0.05) %>% arrange(classicFisher) %>% mutate(ontology = "MF")


# CC ------------------------------------------------------------------------------------------

ob_topGO_selected_ks_cc <- get_gsea_results(ob_selected_genes_pvals, ob_gene2dm_GO_list, mode = "both",
                                            outputTopGOobject = FALSE, nodeSize = 10, ontology = "CC")
# ob_topGO_selected_ks_depleted <- ob_topGO_selected_ks %>% filter(Significant < Expected) %>% 
#    filter(classicFisher < 0.05) %>% arrange(classicFisher)
ob_topGO_selected_ks_cc_enriched <- ob_topGO_selected_ks_cc %>% filter(Significant > Expected) %>% 
   filter(classicFisher < 0.05) %>% arrange(classicFisher) %>% mutate(ontology = "CC")


# deg_genes <- genes$Osmia_bicornis$pvalue
# names(deg_genes) <- genes$Osmia_bicornis$geneID
# deg_genes <- na.omit(deg_genes)
# 
# ob_topGO_sexes <- get_gsea_results(deg_genes, ob_gene2dm_GO_list, mode = "both",
#                                    outputTopGOobject = FALSE, nodeSize = 10)

ob_topGO_merged <- bind_rows(ob_topGO_selected_ks_bp_enriched[1:30, ], 
                             ob_topGO_selected_ks_mf_enriched[1:4, ],
                             ob_topGO_selected_ks_cc_enriched[1:4, ]) %>% 
   mutate(stat = -log2(as.numeric(classicFisher))) %>%
   mutate(Term = paste0(ifelse(nchar(Term) < 67, Term, paste0(str_sub(Term, end = 67), "[...]")), " (", Significant, ")")) %>%
   dplyr::filter(!duplicated(Term)) %>% 
   mutate(TermAsFactor = factor(Term, levels = rev(Term)))

pos_selection_go_plot <- ggplot(ob_topGO_merged, aes(x = TermAsFactor, y = stat, fill = ontology)) + geom_col() + 
   xlab("GO term") + ylab("-log2(p)") + 
   geom_hline(yintercept = c(-log(0.05, 2)), linetype = "dashed") + coord_flip() #+ 
#   scale_y_discrete(limits = rev(levels(TermAsFactor)))

ggsave("04_results/PRJNA285788/pos_selection_go_plot.png", pos_selection_go_plot)



# run again with additional filter LFC > 0 and LFC < 0 (conservation of female-biased and male-biased genes)

shared_among_all_bees <- single_copy_ogs
for (gene_df in genes) {
   selected <- gene_df %>% filter(episDivSelection) %>% drop_na() %>% pull(OG)
   shared_among_all_bees <- intersect(shared_among_all_bees, selected)
}

am_ogs_diff_expr <- genes$Apis_mellifera %>% filter(padj < 0.05) %>% drop_na() %>% pull(OG) 
ob_ogs_diff_expr <- genes$Osmia_bicornis %>% filter(padj < 0.05) %>% drop_na() %>% pull(OG)
bt_ogs_diff_expr <- genes$Bombus_terrestris %>% filter(padj < 0.05) %>% drop_na() %>% pull(OG)
shared_ogs_diff_expr <- intersect(intersect(am_ogs_diff_expr, ob_ogs_diff_expr), bt_ogs_diff_expr)

am_selected_ogs <- genes$Apis_mellifera %>% filter(episDivSelection) %>% drop_na() %>% pull(OG)
ob_selected_ogs <- genes$Osmia_bicornis %>% filter(episDivSelection) %>% drop_na() %>% pull(OG)
bt_selected_ogs <- genes$Bombus_terrestris %>% filter(episDivSelection) %>% drop_na() %>% pull(OG)
ol_selected_ogs <- genes$Osmia_lignaria %>% filter(episDivSelection) %>% drop_na() %>% pull(OG)
shared_selected_ogs <- intersect(
   intersect(intersect(am_selected_ogs, ob_selected_ogs), bt_selected_ogs), ol_selected_ogs)

am_selected_ogs_diff_expr <- genes$Apis_mellifera %>% filter(padj < 0.05, episDivSelection) %>% drop_na() %>% pull(OG)
ob_selected_ogs_diff_expr <- genes$Osmia_bicornis %>% filter(padj < 0.05, episDivSelection) %>% drop_na() %>% pull(OG)
bt_selected_ogs_diff_expr <- genes$Bombus_terrestris %>% filter(padj < 0.05, episDivSelection) %>% drop_na() %>% pull(OG)
shared_selected_ogs_diff_expr <- intersect(
   intersect(am_selected_ogs_diff_expr, ob_selected_ogs_diff_expr), bt_selected_ogs_diff_expr)

am_ogs_diff_expr_low_div <- genes$Apis_mellifera %>% filter(lowDivGene, padj < 0.05) %>% drop_na() %>% pull(OG)
bt_ogs_diff_expr_low_div <- genes$Bombus_terrestris %>% filter(lowDivGene, padj < 0.05) %>% drop_na() %>% pull(OG)
ob_ogs_diff_expr_low_div <- genes$Osmia_bicornis %>% filter(lowDivGene, padj < 0.05) %>% drop_na() %>% pull(OG)
shared_ogs_diff_expr_low_div <- intersect(
   intersect(am_ogs_diff_expr_low_div, bt_ogs_diff_expr_low_div), ob_ogs_diff_expr_low_div)
