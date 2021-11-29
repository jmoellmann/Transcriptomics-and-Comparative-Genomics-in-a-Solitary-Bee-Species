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

in_dir <- paste0(project_dir, "/04_results/PRJNA533756/after_trimming/quantification/STAR")
out_dir <- paste0(project_dir, "/04_results/PRJNA533756/after_trimming/diff_expr/STAR")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(paste0(in_dir, "/merged"))

# get metadata --------------------------------------------------------------------------------

metadata <- read_csv(paste0(project_dir, "/01_data/01_rnaseq_datasets/PRJNA533756_warner/accession/SraRunTable.csv"))

experim_groups <- c("adult male", "mated queen", "virgin queen", "worker forager", "worker nurse")

metadata <- metadata %>% dplyr::select(Run, SAMPLE_TYPE, Replicate, Tissue) %>% 
  dplyr::rename(Group = SAMPLE_TYPE) %>% filter(Group %in% experim_groups) %>% 
  mutate(Group = str_replace(Group, " ", "_")) %>% 
  mutate(Sex = ifelse(Group == "adult_male", "male", "female")) %>% 
  mutate(Name = paste0(Run, "_", Group, "_", Replicate, "_", Tissue)) %>% mutate_all(as.factor) 

experim_groups <- str_replace(experim_groups, " ", "_")

# Merge tissue read counts

files <- file.path(in_dir, paste0(metadata$Run, "_ReadsPerGene.out.tab"))

collapsed_metadata <- metadata %>% group_by(Group, Replicate) %>% 
  summarise(Group, Replicate, Run = paste(as.character(Run), collapse = "_")) %>% 
  ungroup() %>% 
  mutate(Name = paste0(Group, "_", Replicate)) %>% 
  mutate(Sex = ifelse(Group == "adult_male", "male", "female")) %>% mutate_all(as.factor) %>% 
  distinct()
  

hist_plots <- list()
for(runs in collapsed_metadata %>% pull(Runs)){
  runs_list <- str_split(runs, "_")[[1]]
  count_matrices <- list()
  for(run in runs_list){
    tmp_df <- read_tsv(paste0(
      in_dir, "/", run, "_ReadsPerGene.out.tab"), skip = 4, col_names = FALSE, col_types = "cddd")
    tmp_matrix <- as.matrix(tmp_df[, c(2,3,4)])
    row.names(tmp_matrix) <- as.vector(unlist(tmp_df[, 1]))
    count_matrices[[run]] <- tmp_matrix
  }
  new_matrix <- count_matrices[[1]] + count_matrices[[2]] + count_matrices[[3]]
  hist_plots[[runs]] <- ggplot(as_tibble(new_matrix), aes(x = X2)) + geom_density() + xlim(0,3000) + ylim(0,0.0025)
  write.table(x = new_matrix, file = paste0(in_dir, "/merged/", runs, "_ReadsPerGene.out.tab"), 
            sep = "\t", quote = FALSE, col.names = FALSE, row.names = TRUE)
}

# DEA, pairwise (Wald test) -------------------------------------------------------------------

dds_wald <- list()
for(tissue in c("Abdomen", "Thorax", "Head")){
  metadata_tmp <- metadata %>% filter(Tissue == tissue) 
  
  dds_wald[[tissue]] <- de_analysis(in_dir, metadata_tmp, design = ~ Group)
}

# 
# # DEA, pairwise (Wald test), no mated queens --------------------------------------------------
# 
# dds_wald_no_mq <- list()
# dds_stats_no_mq <- list()
# for(tissue in c("Abdomen", "Thorax", "Head")){
#   metadata_tmp <- metadata %>% filter(Tissue == tissue, Group != "mated_queen") 
#   
#   experim_groups_no_mq <- experim_groups[experim_groups != "mated_queen"]
#   
#   dds_wald_no_mq[[tissue]] <- de_analysis(in_dir, metadata_tmp, tx2gene, design = ~ Group)
#   dds_stats_no_mq[[tissue]] <- deg_stats(experim_groups_no_mq, dds_wald_no_mq[[tissue]], alpha = 0.05)
# }
# 
# 
# # DEA, all groups (LRT test) ------------------------------------------------------------------
# 
# dds_lrt <- list()
# for(tissue in c("Abdomen", "Thorax", "Head")){
#   metadata_tmp <- metadata %>% filter(Tissue == tissue) 
#   
#   dds_lrt[[tissue]] <- de_analysis(in_dir, metadata_tmp, tx2gene, design = ~ Group, 
#                                    test = "LRT", full = ~ Group, reduced = ~ 1)
# }
# 
# 
# # DEA, all groups (LRT test), no mated queens -------------------------------------------------
# 
# dds_lrt_no_mq <- list()
# for(tissue in c("Abdomen", "Thorax", "Head")){
#   metadata_tmp <- metadata %>% filter(Tissue == tissue, Group != "mated_queen") 
#   
#   dds_lrt_no_mq[[tissue]] <- de_analysis(in_dir, metadata_tmp, tx2gene, design = ~ Group, 
#                                    test = "LRT", full = ~ Group, reduced = ~ 1)
# }


# DEA, sex ------------------------------------------------------------------------------------

# dds_sex <- list()
# for(tissue in c("Abdomen", "Thorax", "Head")){
#   metadata_tmp <- metadata %>% filter(Tissue == tissue) 
#   
#   dds_sex[[tissue]] <- de_analysis(in_dir, metadata_tmp, design = ~ Sex, 
#                                    test = "LRT", full = ~ Sex, reduced = ~ 1)
# }


# DEA, sex, no mated queens -------------------------------------------------------------------

dds_sex_no_mq <- list()
for(tissue in c("Abdomen", "Thorax", "Head")){
  metadata_tmp <- metadata %>% filter(Tissue == tissue, Group != "mated_queen") 
  
  dds_sex_no_mq[[tissue]] <- de_analysis(in_dir, metadata_tmp, design = ~ Sex, 
                                         test = "LRT", full = ~ Sex, reduced = ~ 1)
}


# DEA, caste, virgin queens -------------------------------------------------------------------

dds_caste_virgin <- list()
for(tissue in c("Abdomen", "Thorax", "Head")){
  allowed_groups <- c("worker_forager", "worker_nurse", "virgin_queen")
  metadata_tmp <- metadata %>% filter(Tissue == tissue, Group %in% allowed_groups) %>% 
    mutate(Caste = ifelse(Group == "virgin_queen", "queen", "worker")) %>% mutate_all(as.factor)
  
  dds_caste_virgin[[tissue]] <- de_analysis(in_dir, metadata_tmp, design = ~ Caste, 
                                            test = "LRT", full = ~ Caste, reduced = ~ 1)
}



# DEA, caste, mated queens --------------------------------------------------------------------

# dds_caste_mated <- list()
# for(tissue in c("Abdomen", "Thorax", "Head")){
#   allowed_groups <- c("worker_forager", "worker_nurse", "mated_queen")
#   metadata_tmp <- metadata %>% filter(Tissue == tissue, Group %in% allowed_groups) %>% 
#     mutate(Caste = ifelse(Group == "mated_queen", "queen", "worker")) %>% mutate_all(as.factor)
#   
#   dds_caste_mated[[tissue]] <- de_analysis(in_dir, metadata_tmp, design = ~ Caste, 
#                                            test = "LRT", full = ~ Caste, reduced = ~ 1)
# }


# DEA, tissue ---------------------------------------------------------------------------------

dds_tissue <- de_analysis(in_dir, metadata, design = ~ Tissue)

# Whole body (merged tissues) sex analysis

tmp_collapsed_metadata <- collapsed_metadata %>% filter(Group != "mated_queen")
dds_merged_sex <- de_analysis(paste0(in_dir, "/merged"), tmp_collapsed_metadata, design = ~ Sex, test = "LRT", 
            full = ~ Sex, nskip = 0, reduced = ~ 1)

# save ----------------------------------------------------------------------------------------


saveRDS(dds_wald, file = paste0(out_dir, "/dds_wald.RDS"))
#saveRDS(dds_sex, file = paste0(out_dir, "/dds_sex.RDS"))
saveRDS(dds_sex_no_mq, file = paste0(out_dir, "/dds_sex_no_mq.RDS"))
#saveRDS(dds_caste_mated, file = paste0(out_dir, "/dds_caste_mated.RDS"))
saveRDS(dds_caste_virgin, file = paste0(out_dir, "/dds_caste_virgin.RDS"))
saveRDS(dds_tissue, file = paste0(out_dir, "/dds_tissue.RDS"))
saveRDS(dds_merged_sex, file = paste0(out_dir, "/dds_merged_sex.RDS"))

#saveRDS(dds_wald_no_mq, file = paste0(out_dir, "/dds_wald_no_mq.RDS"))
#saveRDS(dds_lrt, file = paste0(out_dir, "/dds_lrt.RDS"))
#saveRDS(dds_lrt_no_mq, file = paste0(out_dir, "/dds_lrt_no_mq.RDS"))