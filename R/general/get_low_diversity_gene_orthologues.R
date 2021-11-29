library(tidyverse)

project_dir <- "/media/jannik/Samsung_T5/masters_project/01_main_analysis/"
setwd(project_dir)

get_low_div_genes <- function(species_name){
   abbrev <- paste0(str_sub(str_to_lower(str_split(species_name, "_")[[1]][1]), 1, 1), 
                    str_sub(str_split(species_name, "_")[[1]][2], 1, 1))
   bt_protein2species_protein <- orthofinder_res_bees %>% 
      dplyr::select(OG, Bombus_terrestris, species_name) %>%
      tidyr::separate_rows(species_name, sep = ", ") %>%
      tidyr::separate_rows(Bombus_terrestris, sep = ", ") %>%
      dplyr::rename("Bt_ProteinID" = "Bombus_terrestris", "Species_ProteinID" = species_name) %>%
      dplyr::select("Bt_ProteinID", "Species_ProteinID") %>% distinct() %>% drop_na()
   
   species_gene2protein <- read_tsv(paste0("01_data/03_reference/annotation/", abbrev, "_gene2protein.tsv"),
                               col_names = c("Species_GeneID", "Species_ProteinID"))
   bt_gene2species_gene <- bt_protein2species_protein %>% left_join(bt_gene2protein) %>% 
      left_join(species_gene2protein) %>% dplyr::select(Bt_GeneID, Species_GeneID) %>% drop_na()
   low_div_genes_orthofinder <- bt_gene2species_gene %>% filter(Bt_GeneID %in% bt_low_div_genes) %>% 
      pull(Species_GeneID)
   return(low_div_genes_orthofinder)
}

#insect_orthologues_dir <- 
#   "results/PRJNA285788/immune_gene_identification/orthofinder_insects/Phylogenetic_Hierarchical_Orthogroups"
bee_orthologues_dir <- 
   "04_results/PRJNA285788/immune_gene_identification/orthofinder_bees/Phylogenetic_Hierarchical_Orthogroups"

orthofinder_res_bees <- readr::read_tsv(paste0(bee_orthologues_dir, "/N0.tsv")) 

ob_protein2bt_protein <- orthofinder_res_bees %>% 
   dplyr::select(OG, Bombus_terrestris, Osmia_bicornis) %>%
   tidyr::separate_rows(Osmia_bicornis, sep = ", ") %>% 
   tidyr::separate_rows(Bombus_terrestris, sep = ", ") %>% 
   dplyr::rename("Ob_ProteinID" = "Osmia_bicornis", "Bt_ProteinID" = "Bombus_terrestris") %>% 
   dplyr::select("Ob_ProteinID", "Bt_ProteinID") %>% distinct() %>% drop_na()

ob_protein2am_protein <- orthofinder_res_bees %>% 
   dplyr::select(OG, Apis_mellifera, Osmia_bicornis) %>%
   tidyr::separate_rows(Osmia_bicornis, sep = ", ") %>%
   tidyr::separate_rows(Apis_mellifera, sep = ", ") %>%
   dplyr::rename("Ob_ProteinID" = "Osmia_bicornis", "Am_ProteinID" = "Apis_mellifera") %>%
   dplyr::select("Ob_ProteinID", "Am_ProteinID") %>% distinct() %>% drop_na()

bt_protein2am_protein <- orthofinder_res_bees %>% 
   dplyr::select(OG, Bombus_terrestris, Apis_mellifera) %>%
   tidyr::separate_rows(Apis_mellifera, sep = ", ") %>%
   tidyr::separate_rows(Bombus_terrestris, sep = ", ") %>%
   dplyr::rename("Bt_ProteinID" = "Bombus_terrestris", "Am_ProteinID" = "Apis_mellifera") %>%
   dplyr::select("Bt_ProteinID", "Am_ProteinID") %>% distinct() %>% drop_na()

ob_gene2protein <- read_tsv("01_data/03_reference/annotation/ob_gene2protein.tsv", 
                               col_names = c("Ob_GeneID", "Ob_ProteinID"))

bt_gene2protein <- read_tsv("01_data/03_reference/annotation/bt_gene2protein.tsv", 
                               col_names = c("Bt_GeneID", "Bt_ProteinID"))

am_gene2protein <- read_tsv("01_data/03_reference/annotation/am_gene2protein_v2.tsv", 
                               col_names = c("Am_GeneID", "Am_ProteinID"))


bt_gene2am_gene <- bt_protein2am_protein %>% left_join(bt_gene2protein) %>% 
   left_join(am_gene2protein) %>% dplyr::select(Bt_GeneID, Am_GeneID) %>% drop_na()

ob_gene2am_gene <- ob_protein2am_protein %>% left_join(ob_gene2protein) %>% 
   left_join(am_gene2protein) %>% dplyr::select(Ob_GeneID, Am_GeneID) %>% drop_na()

ob_gene2bt_gene <- ob_protein2bt_protein %>% left_join(ob_gene2protein) %>% 
   left_join(bt_gene2protein) %>% dplyr::select(Ob_GeneID, Bt_GeneID) %>% drop_na()


bt_low_div_genes <- read_table("01_data/03_reference/low_diversity_genes/bt_low_diversity_genes.txt", 
                                     col_names = FALSE)[[1]]
am_low_div_genes <- read_table("01_data/03_reference/low_diversity_genes/am_low_diversity_genes.txt",
                                     col_names = FALSE)[[1]]

am_low_div_genes_orthofinder <- bt_gene2am_gene %>% filter(Bt_GeneID %in% bt_low_div_genes) %>% pull(Am_GeneID)
ob_low_div_genes_orthofinder <- ob_gene2bt_gene %>% filter(Bt_GeneID %in% bt_low_div_genes) %>% pull(Ob_GeneID)


# Other species -------------------------------------------------------------------------------

low_div_genes <- list()
species_names <- c("Megalopta_genalis", "Dufourea_novaeangliae", "Ceratina_calcarata", 
                   "Eufriesea_mexicana", "Habropoda_laboriosa", "Megachile_rotundata", 
                   "Nomia_melanderi", "Osmia_lignaria")
for(species_name in species_names){
   low_div_genes[[species_name]] <- get_low_div_genes(species_name)
   abbrev <- paste0(str_sub(str_to_lower(str_split(species_name, "_")[[1]][1]), 1, 1),
                                 str_sub(str_split(species_name, "_")[[1]][2], 1, 1))
   write.table(low_div_genes[[species_name]], file = paste0("01_data/03_reference/low_diversity_genes/",
                                                        abbrev, "_low_diversity_genes_orthofinder.txt"),
               quote = FALSE, sep = "", row.names = FALSE, col.names = FALSE)

}

write.table(am_low_div_genes_orthofinder, file = "01_data/03_reference/low_diversity_genes/am_low_diversity_genes_orthofinder.txt",
            quote = FALSE, sep = "", row.names = FALSE, col.names = FALSE)
write.table(ob_low_div_genes_orthofinder, file = "01_data/03_reference/low_diversity_genes/ob_low_diversity_genes_orthofinder.txt",
            quote = FALSE, sep = "", row.names = FALSE, col.names = FALSE)
