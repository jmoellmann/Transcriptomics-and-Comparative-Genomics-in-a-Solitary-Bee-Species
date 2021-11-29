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

in_dir <- paste0(project_dir, "/results/PRJNA285788/pesticide_exposure/quantification/STAR")
out_dir <- paste0(project_dir, "/results/PRJNA285788/pesticide_exposure/diff_expr/STAR")

dir.create(out_dir, recursive = TRUE)

# get metadata --------------------------------------------------------------------------------

metadata <- read.csv(paste0(project_dir, "/data/datasets/PRJNA285788/accession/SraRunTable.csv"), 
                     header = TRUE)

experim_groups <- c("Control", "Imidacloprid", "Thiacloprid")

metadata <- metadata %>% dplyr::select(Run, LibrarySelection, isolate) %>% 
  dplyr::mutate(Group = unlist(lapply(
    stringr::str_split(isolate, pattern = " "), function(x) x[1]))) %>% 
  dplyr::filter(LibrarySelection == "cDNA") %>% dplyr::select(Run, Group) %>% 
  dplyr::mutate_all(as.factor) 

immune_genes <- read.table("results/PRJNA285788/immune_gene_identification/immune_genes/ob_merged_imm_genes.txt")[[1]]

osmia_specific_genes <- read_tsv("results/PRJNA285788/immune_gene_identification/osmia_specific_genes.tsv") %>% 
  pull(GeneID)

# DEA, (Wald test) -------------------------------------------------------------------

dds_wald <- de_analysis(in_dir, metadata, design = ~ Group, test = "Wald")

res <- list()

for(treatment in c("Thiacloprid", "Imidacloprid")){
  res[[treatment]] <- results(dds_wald, contrast = c("Group", treatment, "Control")) %>% 
    tibble::as_tibble() %>% 
    dplyr::mutate(geneID = row.names(results(dds_wald, contrast = c("Group", treatment, "Control")))) %>% 
    dplyr::arrange(geneID)
  
  geneID2geneDescription <- read_delim("data/reference/annotation/ob_geneID2geneDescription.txt", delim = ";", 
                                       col_names = c("geneID", "geneDescription")) %>% 
    mutate(geneID = str_extract(geneID, "LOC[0-9]+"), 
           geneDescription = str_remove(geneDescription, "product"),
           geneDescription = str_remove_all(geneDescription, "\""),
           geneDescription = str_remove(geneDescription, "\\s+")) %>% 
    dplyr::filter(!is.na(geneID), !is.na(geneDescription)) %>% 
    group_by(geneID) %>% summarise(geneDescription = min(geneDescription)) %>% 
    right_join(res[[treatment]]) %>% arrange(geneID)
  
  res[[treatment]] <- res[[treatment]] %>% mutate(geneDescription = geneID2geneDescription$geneDescription) %>% 
    dplyr::select(geneID, geneDescription, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj) %>% 
    arrange(padj)
  
  write_csv(res[[treatment]], paste0("results/PRJNA285788/pesticide_exposure/diff_expr/STAR/DESeq2_results_", 
                        tolower(treatment), ".csv"))
  
  write_csv(res[[treatment]] %>% dplyr::filter(geneID %in% immune_genes), 
            paste0("results/PRJNA285788/pesticide_exposure/diff_expr/STAR/DESeq2_results_immune_genes_", 
                   tolower(treatment), ".csv"))
  
  write_csv(res[[treatment]] %>% dplyr::filter(geneID %in% osmia_specific_genes), 
            paste0("results/PRJNA285788/pesticide_exposure/diff_expr/STAR/DESeq2_results_osmia_specific_genes_", 
            tolower(treatment), ".csv"))
}

# save ----------------------------------------------------------------------------------------

saveRDS(dds_wald, file = paste0(out_dir, "/dds_wald.RDS"))
