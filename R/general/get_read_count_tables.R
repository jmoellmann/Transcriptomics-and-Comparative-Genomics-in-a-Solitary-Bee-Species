library(tximport)
library(dplyr)
library(readr)

# start of analysis ---------------------------------------------------------------------------

project_dir <- "/media/jannik/Samsung_T5/masters_project/ApisTranscriptomics"
setwd(project_dir)

# get metadata --------------------------------------------------------------------------------

metadata <- read.csv(paste0(project_dir, "/data/datasets/PRJNA285788/accession/SraRunTable.csv"), 
                     header = TRUE)

metadata <- metadata[metadata$Assay.Type == "RNA-Seq", ] %>% dplyr::select("Run", "Assay.Type") %>% 
  dplyr::mutate_all(as.factor) 

# STAR ----------------------------------------------------------------------------------------

sex_differences_star <- tibble::as_tibble(counts(readRDS(
  "results/PRJNA285788/sex_differences/diff_expr/STAR/dds_lrt.RDS")), rownames = "GeneID")
pesticide_exposure_star <- tibble::as_tibble(counts(readRDS(
  "results/PRJNA285788/pesticide_exposure/diff_expr/STAR/dds_wald.RDS")), 
  rownames = "GeneID")

read_counts_star <- full_join(sex_differences_star, pesticide_exposure_star) %>% 
  relocate("GeneID", sort(metadata$Run)) %>% arrange(GeneID)

write_csv(read_counts_star, "data/datasets/PRJNA285788/read_counts_star.csv")

# Kallisto ------------------------------------------------------------------------------------

sex_differences_kallisto <- tibble::as_tibble(counts(readRDS(
  "results/PRJNA285788/sex_differences/diff_expr/kallisto/dds_lrt.RDS")), rownames = "GeneID")
pesticide_exposure_kallisto <- tibble::as_tibble(counts(readRDS(
  "results/PRJNA285788/pesticide_exposure/diff_expr/kallisto/dds_wald.RDS")), 
  rownames = "GeneID")

read_counts_kallisto <- full_join(sex_differences_kallisto, pesticide_exposure_kallisto) %>% 
  relocate("GeneID", sort(metadata$Run)) %>% arrange(GeneID)

write_csv(read_counts_kallisto, "data/datasets/PRJNA285788/read_counts_kallisto.csv")
