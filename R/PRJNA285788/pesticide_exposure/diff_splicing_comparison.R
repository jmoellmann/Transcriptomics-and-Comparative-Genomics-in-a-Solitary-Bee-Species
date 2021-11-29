library(tidyverse)

project_dir <- "/media/jannik/Samsung_T5/masters_project/ApisTranscriptomics/"
setwd(project_dir)

in_dir <- paste0(project_dir, "/results/PRJNA285788/pesticide_exposure/diff_splicing/")

# get metadata --------------------------------------------------------------------------------

metadata <- read.csv(paste0(project_dir, "/data/datasets/PRJNA285788/accession/SraRunTable.csv"), 
                     header = TRUE)

metadata <- metadata %>% dplyr::select(Run, LibrarySelection, sex) %>% 
  dplyr::rename(Sex = sex) %>% dplyr::filter(LibrarySelection == "PolyA") %>% dplyr::mutate_all(as.factor) 

immune_genes <- read.table(
  "results/PRJNA285788/immune_gene_identification/immune_genes/ob_merged_imm_genes.txt")[["V1"]]

all_genes <- row.names(DESeq2::results(readRDS(paste0(
  "results/PRJNA285788/sex_differences/diff_expr/STAR/dds_lrt.RDS"))))

# splicing ------------------------------------------------------------------------------------

res <- list()

geneID2geneDescription <- read_delim("data/reference/annotation/ob_geneID2geneDescription.txt", delim = ";", 
                                     col_names = c("geneID", "geneDescription")) %>% 
  mutate(geneID = str_extract(geneID, "LOC[0-9]+"), 
         geneDescription = str_remove(geneDescription, "product"),
         geneDescription = str_remove_all(geneDescription, "\""),
         geneDescription = str_remove(geneDescription, "\\s+")) %>% 
  dplyr::filter(!is.na(geneID), !is.na(geneDescription)) %>% 
  group_by(geneID) %>% summarise(geneDescription = min(geneDescription)) %>% 
  dplyr::rename(GeneID = geneID)

for(treatment in c("imidacloprid", "thiacloprid")){
  vect_binom.test <- Vectorize(binom.test)
  
  summary <- read_tsv(paste0(in_dir, treatment, "/rMATS_results/summary.txt"))
  
  new_col <- c()
  for(row in 1:dim(summary)[1]){
    new_col <- c(new_col, binom.test(x = pull(summary[row, "SigEventsJCSample1HigherInclusion"]), 
                                     n = pull(summary[row, "SignificantEventsJC"], p = 0.5))$p.value)
  }
  summary$pvals <- new_col
  summary$padj <- p.adjust(summary$pvals, method = "BH")
  
  write_csv(summary, paste0("results/PRJNA285788/pesticide_exposure/diff_splicing/", treatment, 
                            "/summary.csv"))
  
  res_temp <- tibble()
  for(event in c("A3SS", "A5SS", "MXE", "RI", "SE")){
    df <- read_tsv(paste0("results/PRJNA285788/pesticide_exposure/diff_splicing/", tolower(treatment), 
                          "/rMATS_results/", event, ".MATS.JC.txt")) %>% mutate(event = event)
    res_temp <- res_temp %>% bind_rows(df)
  }
  
  
  #res %>% filter(FDR < 0.05) %>% group_by(event) %>% summarise(n = length(PValue))
  
  res_all <- res_temp %>% group_by(GeneID) %>% summarise(minPval = min(PValue), minFDR = min(FDR))
  
  res_temp <- res_temp %>% left_join(geneID2geneDescription) %>% arrange(FDR) %>% 
    dplyr::select(GeneID, geneDescription, chr, strand, event, PValue:IncLevelDifference, 
                  longExonStart_0base:flankingEE, IJC_SAMPLE_1:SkipFormLen,
                  `1stExonStart_0base`:downstreamEE, 
                  riExonStart_0base:exonEnd)
  
  write_csv(res_all, paste0("results/PRJNA285788/pesticide_exposure/diff_splicing/", 
                            tolower(treatment), "res_all.csv"))
  write_csv(res_temp, paste0("results/PRJNA285788/pesticide_exposure/diff_splicing/", 
                             tolower(treatment), "rMATS_results.csv"))
  
  res[[treatment]] <- res_all
  
  res_dsg <- res_all %>% filter(minFDR < 0.05)
  res_immune_dsg <- res_dsg %>% filter(GeneID %in% immune_genes)
  
  test_matrix <- matrix(data = c(length(immune_genes), length(res_immune_dsg$GeneID), 
                                 length(all_genes), length(res_dsg$GeneID)), 
                        nrow = 2, ncol = 2)
  
  f.test <- fisher.test(test_matrix)
  
  out_file <- paste0("results/PRJNA285788/pesticide_exposure/diff_splicing/", treatment, "/stats.txt")
  sink(out_file, append = FALSE)
  cat(paste("Genes:", test_matrix[1,2], "\n"))
  cat(paste("Immune Genes:", test_matrix[1,1], "\n"))
  cat(paste("DSGs:", test_matrix[2,2], "\n"))
  cat(paste("Immune DSGs:", test_matrix[2,1], "\n"))
  cat("-------------------------------------------------", "\n")
  cat(paste0("Share of immune genes diff. spliced: ", 
             round((test_matrix[2,1] / test_matrix[1,1]), 4) * 100, "%", "\n"))
  cat(paste0("Share of all genes diff. spliced.: ", 
             round((test_matrix[2,2] / test_matrix[1,2]), 4) * 100, "%", "\n"))
  cat("-------------------------------------------------", "\n")
  cat("\n")
  cat(paste0("Fisher-test p-value: ", round(f.test$p.value, 3), "\n"))
  sink()
  
}
