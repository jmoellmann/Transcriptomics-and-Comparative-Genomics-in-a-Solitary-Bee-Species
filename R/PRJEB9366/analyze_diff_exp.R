library(tidyverse)
library(DESeq2)

de_analysis <- function(in_dir, metadata, design, test = "Wald", full = NULL, reduced = ~ 1, 
                        nskip = 4){
   
   files <- file.path(in_dir, paste0(metadata$Run, "_ReadsPerGene.out.tab"))
   count_files <- lapply(files, read_tsv, skip = nskip, col_names = FALSE, col_types = "cddd")
   
   count_matrix <- as.data.frame(sapply(count_files, function(x) x[, 2]))
   colnames(count_matrix) <- metadata$Name
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

deg_stats <- function(experim_groups, dds_object, alpha){
   n = length(experim_groups)
   m <- matrix(nrow = n, ncol = n)
   for(i in 1:(n - 1)){
      for(j in (i + 1):n){
         res <- results(dds_object, contrast = c("Group", experim_groups[i], experim_groups[j]))
         padj_nona <- res$padj[!is.na(res$padj)]
         n_deg = length(padj_nona[padj_nona < alpha])
         
         m[i,j] <- n_deg
      }
   }
   colnames(m) <- experim_groups
   rownames(m) <- experim_groups
   return(m)
}


# start of analysis ---------------------------------------------------------------------------

project_dir <- "/media/jannik/Samsung_T5/masters_project/01_main_analysis/"

in_dir <- paste0(project_dir, "/04_results/PRJEB9366/quantification/STAR")
out_dir <- paste0(project_dir, "/04_results/PRJEB9366/diff_expr/STAR")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(paste0(in_dir, "/merged"))

# get metadata --------------------------------------------------------------------------------

metadata <- read_tsv(paste0(project_dir, "/01_data/01_rnaseq_datasets/PRJEB9366/accession/ena_accession.tsv"))

experim_groups <- c("MA", "WAR", "WAU", "VQ")

metadata <- metadata %>% dplyr::select(run_accession, submitted_ftp, group) %>% 
   filter(group %in% experim_groups) %>% 
   mutate(Sex = ifelse(group == "MA", "male", "female")) %>% 
   mutate(Run = sapply(str_split(sapply(str_split(submitted_ftp, "/"), function(x) x[6]), "\\."), 
                       function(x) x[1])) %>% 
   mutate(Name = Run) %>% mutate_all(as.factor) 

# DEA, pairwise (Wald test) -------------------------------------------------------------------

dds_wald <- de_analysis(in_dir, metadata, design = ~ group)

# DEA, LRT ------------------------------------------------------------------------------------

dds_lrt <- de_analysis(in_dir, metadata, design = ~ Sex, test = "LRT", full = ~ Sex)

# save ----------------------------------------------------------------------------------------

saveRDS(dds_wald, file = paste0(out_dir, "/dds_wald.RDS"))
saveRDS(dds_lrt, file = paste0(out_dir, "/dds_lrt.RDS"))
