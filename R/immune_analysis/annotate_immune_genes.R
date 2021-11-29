library(tidyverse)

get_type <- function(y, dm_unique_genes, bt_unique_genes, dm_bt_overlap_genes){
   types <- c()
   for(x in y){
      if(x %in% dm_unique_genes){
         types <- c(types, "dm-unique")
      } else if(x %in% bt_unique_genes){
         types <- c(types, "bt-unique")
      } else if(x %in% dm_bt_overlap_genes){
         types <- c(types, "dm-bt-shared")
      } else{
         types <- c(types, NA)
      }
   }
   return(types)
}

get_parent_group <- function(table){
   for(i in 1:nrow(table)){
      j = i
      while(as.vector(!is.na(table[j, "ParentGroupID"]))){
         j = which(table[["GroupID"]] == table[j, "ParentGroupID"])[1]
      }
      table[i, "Group"] = table[j, "Group"]
   }
   return(table)
}

project_dir <- "/media/jannik/Samsung_T5/masters_project/ApisTranscriptomics/"
setwd(project_dir)

immune_genes <- read.table("results/PRJNA285788/immune_gene_identification/immune_genes/ob_merged_imm_genes.txt")
colnames(immune_genes) <- c("geneID")

dm_genes <- read.table("results/PRJNA285788/immune_gene_identification/immune_genes/ob_dm_imm_genes.txt")[[1]]
bt_genes <- read.table("results/PRJNA285788/immune_gene_identification/immune_genes/ob_bt_imm_genes.txt")[[1]]

dm_bt_overlap_genes <- intersect(dm_genes, bt_genes)
dm_unique_genes <- dm_genes[!(dm_genes %in% dm_bt_overlap_genes)]
bt_unique_genes <- bt_genes[!(bt_genes %in% dm_bt_overlap_genes)]

dm_gene_group2dm_gene <- read_tsv("data/reference/annotation/gene_group_data_fb_2021_03.tsv", 
                                  skip = 8) %>% 
   dplyr::select('## FB_group_id', FB_group_name, Parent_FB_group_id, Group_member_FB_gene_id) %>% 
   dplyr::rename(Dm_GeneID = Group_member_FB_gene_id, Group = FB_group_name,
                 GroupID = '## FB_group_id', ParentGroupID = Parent_FB_group_id)

ob_gene2dm_gene <- read_csv("results/PRJNA285788/immune_gene_identification/ob_gene2dm_gene.csv")

ob_gene2gene_group <- get_parent_group(as.data.frame(dm_gene_group2dm_gene)) %>% as_tibble() %>% 
   dplyr::select(Group, Dm_GeneID) %>% left_join(ob_gene2dm_gene) %>% dplyr::select(Group, Ob_GeneID)

ob_gene2chr_start_end <- read_tsv("data/reference/annotation/ob_gene2chr_start_end.tsv", 
                                  col_names = c("geneID", "chromosome", "strand", "start", "end")) %>% 
   mutate(chromosome = ifelse(str_detect(chromosome, "CAJRAL"), "unplaced", chromosome))

immune_genes_annotated <- read_tsv("data/reference/annotation/ob_gene2gene_description.tsv", 
                                     col_names = c("geneID", "geneDescription")) %>% 
   dplyr::filter(!is.na(geneID), !is.na(geneDescription)) %>% 
   group_by(geneID) %>% summarise(geneDescription = min(geneDescription)) %>% 
   right_join(immune_genes) %>% 
   mutate(origin = get_type(geneID, dm_unique_genes, bt_unique_genes, dm_bt_overlap_genes)) %>% 
   left_join(ob_gene2gene_group, by = c("geneID" = "Ob_GeneID")) %>% group_by(geneID) %>% 
   summarise(geneID, geneDescription, origin, group = paste(str_to_title(unique(Group)), 
                                                            collapse = ", ")) %>% 
   distinct() %>% arrange(origin) %>% left_join(ob_gene2chr_start_end)

write_csv(immune_genes_annotated, "results/PRJNA285788/immune_gene_identification/immune_genes_annotated.csv")
