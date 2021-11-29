library(dplyr)
library(readr)
library(stringr)
library(DESeq2)
library(tibble)

de_analysis <- function(in_dir, metadata, design, test = "Wald", full = NULL, reduced = ~ 1){
  
  files <- file.path(in_dir, paste0(metadata$Run, "_ReadsPerGene.out.tab"))
  count_files <- lapply(files, read_tsv, skip = 4, col_names = FALSE, col_types = "cddd")
  
  count_matrix <- as.data.frame(sapply(count_files, function(x) x[, 2]))
  colnames(count_matrix) <- metadata$Run
  rownames(count_matrix) <- count_files[1][[1]]$X1
  
  dds <- DESeqDataSetFromMatrix(countData = count_matrix, 
                                colData = metadata, 
                                design = design)
  if(test == "Wald"){
    return(DESeq(dds))
  } else if(test == "LRT"){
    return(DESeq(dds, test = test, full = full, reduced = reduced))
  }
}

# start of analysis ---------------------------------------------------------------------------

project_dir <- "/media/jannik/Samsung_T5/masters_project/ApisTranscriptomics"

setwd(project_dir)

in_dir <- paste0(project_dir, "/results/PRJNA285788/sex_differences/quantification/STAR")
out_dir <- paste0(project_dir, "/results/PRJNA285788/sex_differences/diff_expr/STAR")

dir.create(out_dir, recursive = TRUE)


# get metadata --------------------------------------------------------------------------------

metadata <- read.csv(paste0(project_dir, "/data/datasets/PRJNA285788/accession/SraRunTable.csv"), 
                     header = TRUE)

experim_groups <- c("female", "male")

metadata <- metadata %>% dplyr::select(Run, LibrarySelection, sex) %>% 
  dplyr::rename(Sex = sex) %>% dplyr::filter(LibrarySelection == "PolyA") %>% dplyr::mutate_all(as.factor) 

immune_genes <- read.table("results/PRJNA285788/immune_gene_identification/immune_genes/ob_merged_imm_genes.txt")[[1]]

osmia_specific_genes <- read_tsv("results/PRJNA285788/immune_gene_identification/osmia_specific_genes.tsv") %>% 
  pull(GeneID)

# DEA, (LRT test) -------------------------------------------------------------------

dds_lrt <- de_analysis(in_dir, metadata, design = ~ Sex, test = "LRT", full = ~ Sex)

res <- tibble::as_tibble(results(dds_lrt)) %>% mutate(geneID = row.names(results(dds_lrt))) %>% 
  arrange(geneID)

geneID2geneDescription <- read_delim("data/reference/annotation/ob_geneID2geneDescription.txt", delim = ";", 
           col_names = c("geneID", "geneDescription")) %>% 
  mutate(geneID = str_extract(geneID, "LOC[0-9]+"), 
         geneDescription = str_remove(geneDescription, "product"),
         geneDescription = str_remove_all(geneDescription, "\""),
         geneDescription = str_remove(geneDescription, "\\s+")) %>% 
  dplyr::filter(!is.na(geneID), !is.na(geneDescription)) %>% 
  group_by(geneID) %>% summarise(geneDescription = min(geneDescription)) %>% 
  right_join(res) %>% arrange(geneID)

res <- res %>% mutate(geneDescription = geneID2geneDescription$geneDescription) %>% 
  dplyr::select(geneID, geneDescription, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj) %>% 
  arrange(padj)

write_csv(res, "results/PRJNA285788/sex_differences/diff_expr/STAR/DESeq2_results.csv")

write_csv(res %>% dplyr::filter(geneID %in% immune_genes), 
          "results/PRJNA285788/sex_differences/diff_expr/STAR/DESeq2_results_immune_genes.csv")

write_csv(res %>% dplyr::filter(geneID %in% osmia_specific_genes), 
          "results/PRJNA285788/sex_differences/diff_expr/STAR/DESeq2_results_osmia_specific_genes.csv")

# save ----------------------------------------------------------------------------------------

saveRDS(dds_lrt, file = paste0(out_dir, "/dds_lrt.RDS"))
