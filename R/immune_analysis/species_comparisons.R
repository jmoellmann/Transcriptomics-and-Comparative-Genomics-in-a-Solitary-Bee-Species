library(biomaRt)
library(KEGGREST)
library(EnrichmentBrowser)
library(reshape2)
library(RColorBrewer)
library(tidyverse)

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

wilcox.test_ob <- function(col){
  return(wilcox.test(col, gene_counts$Osmia_bicornis)$p.value)
}

count_names <- function(strings){
  out <- c()
  for(string in strings){
    temp <- str_remove(str_split(string, ",")[[1]], " ")
    temp <- temp[!is.na(temp)]
    out <- c(out, length(temp))
  }
  return(out)
}

join_and_count_unique_names <- function(strings){
  out <- c()
  for(string in strings){
    temp <- str_remove(str_split(string, ",")[[1]], " ")
    temp <- temp[!is.na(temp)]
    out <- c(out, temp)
  }
  return(length(unique(out)))
}

vectorized_grepl <- Vectorize(grepl, vectorize.args = "pattern")


project_dir <- "/media/jannik/Samsung_T5/masters_project/ApisTranscriptomics/"
setwd(project_dir)

immune_prots <- read_tsv(
  "results/PRJNA285788/immune_gene_identification/immune_genes/ob_merged_imm_prots.txt", 
  col_names = FALSE) %>% pull()

immune_genes <- read_tsv(
  "results/PRJNA285788/immune_gene_identification/immune_genes/ob_merged_imm_genes.txt", 
  col_names = FALSE) %>% pull()

# ob2dm_gene <- read_csv("results/PRJNA285788/immune_gene_identification/ob_gene2dm_gene.csv") %>% 
#   mutate(Dm_GeneID = str_replace(Dm_GeneID, "[^0-9]+", "")) %>% 
#   mutate(Dm_GeneID = str_replace(Dm_GeneID, "^0+", "")) %>% 
#   filter(!is.na(Dm_GeneID))

dm_mart <- biomaRt::useEnsemblGenomes(biomart = "metazoa_mart", dataset = "dmelanogaster_eg_gene")
dm_gene2dm_protein <- biomaRt::getBM(c("ensembl_gene_id", "refseq_peptide"), 
                                     mart = dm_mart) %>%
  as_tibble() %>% 
  dplyr::filter(refseq_peptide != "") %>% 
  dplyr::rename("Dm_ProteinID" = "refseq_peptide", "Dm_GeneID" = "ensembl_gene_id")

# dm_immune_genes <- ob2dm_gene %>% filter(Ob_GeneID %in% immune_genes) %>% pull(Dm_GeneID) %>% 
#   unique()
# 
# mask <- sapply(dme_pathways, function(x) any(dm_immune_genes %in% x))
# 
# immune_pathways <- dme_pathways[mask]

immune_orthogroups <- read_tsv(
  "results/PRJNA285788/immune_gene_identification/orthofinder_insects/Phylogenetic_Hierarchical_Orthogroups/N0.tsv") %>% 
  dplyr::select(OG, Osmia_bicornis) %>% 
  filter(apply(vectorized_grepl(immune_prots, Osmia_bicornis), MARGIN = 1, 
               function(x) any(x == TRUE))) %>% pull(OG) %>% unique()

write_lines(file = "results/PRJNA285788/immune_gene_identification/immune_genes/imm_ogs.txt", 
           x = immune_orthogroups)

# dm_pathways <- getGenesets(org = "dme", db = "kegg", cache = TRUE, return.type="list")
# 
# dm_pathway2dm_gene <- as_tibble(do.call(cbind, dm_pathways)) %>% 
#   pivot_longer(cols = 1:length(dm_pathways)) %>% distinct() %>% 
#   dplyr::rename(Dm_GeneID = value, Pathway = name) %>% 
#   mutate(Pathway = str_replace_all(Pathway, "_", " ")) %>% 
#   mutate(Pathway = str_remove(Pathway, "dme[0-9]+")) %>% 
#   mutate(Pathway = str_remove(Pathway, " "))

dm_gene_group2dm_gene <- read_tsv("data/reference/annotation/gene_group_data_fb_2021_03.tsv", 
                                  skip = 8) %>% 
  dplyr::select('## FB_group_id', FB_group_name, Parent_FB_group_id, Group_member_FB_gene_id) %>% 
  dplyr::rename(Dm_GeneID = Group_member_FB_gene_id, Group = FB_group_name,
                GroupID = '## FB_group_id', ParentGroupID = Parent_FB_group_id)

# dm_gene_group2dm_gene_v2 <- dm_gene_group2dm_gene %>% dplyr::filter(is.na(ParentGroupID))

dm_gene_group2dm_gene_v2 <- get_parent_group(as.data.frame(dm_gene_group2dm_gene)) %>% as_tibble() %>% 
  dplyr::select(Group, Dm_GeneID)

dm_gene_group2og <- read_tsv(
  "results/PRJNA285788/immune_gene_identification/orthofinder_insects/Phylogenetic_Hierarchical_Orthogroups/N0.tsv") %>% 
  dplyr::select(OG, Drosophila_melanogaster) %>% 
  dplyr::rename(Dm_ProteinID = Drosophila_melanogaster) %>% 
  separate_rows(Dm_ProteinID, sep = ",") %>% 
  mutate(Dm_ProteinID = str_remove(Dm_ProteinID, " ")) %>%
  inner_join(dm_gene2dm_protein) %>% 
  dplyr::select(OG, Dm_GeneID, Dm_ProteinID) %>% inner_join(dm_gene_group2dm_gene_v2) %>% 
  dplyr::rename(Orthogroup = OG)


# dm_pathway2og <- read_tsv(
#   "results/PRJNA285788/immune_gene_identification/orthofinder_insects/Orthogroups/Orthogroups.tsv") %>% 
#   dplyr::select(Orthogroup, Drosophila_melanogaster) %>% 
#   dplyr::rename(Dm_ProteinID = Drosophila_melanogaster) %>% 
#   separate_rows(Dm_ProteinID, sep = ",") %>% 
#   mutate(Dm_ProteinID = str_remove(Dm_ProteinID, " ")) %>%
#   full_join(dm_gene2dm_protein) %>% 
#   dplyr::select(Orthogroup, Dm_GeneID) %>%
#   mutate(Dm_GeneID = str_replace(Dm_GeneID, "[^0-9]+", "")) %>% 
#   mutate(Dm_GeneID = str_replace(Dm_GeneID, "^0+", "")) %>% full_join(dm_pathway2dm_gene) %>% 
#   drop_na()

# 
# gene_counts <- read_tsv(
#   "results/PRJNA285788/immune_gene_identification/orthofinder_insects/Orthogroups/Orthogroups.GeneCount.tsv")
# gene_counts[gene_counts == 0] <- NA
# 
# 
# means <- as_tibble(apply(gene_counts[,2:22], MARGIN = 2, mean, na.rm = TRUE)) %>% 
#   mutate(Species = colnames(gene_counts[,2:22])) %>% dplyr::rename(MeanGeneCount = value) %>% 
#   dplyr::select(Species, MeanGeneCount) %>% 
#   mutate(pVal = apply(gene_counts[,2:22], MARGIN = 2, wilcox.test_ob)) %>% 
#   mutate(pAdj = p.adjust(pVal)) %>% 
#   arrange(desc(MeanGeneCount))
# 
# means[means$Species == "Osmia_bicornis", "pVal"] <- NA
# means[means$Species == "Osmia_bicornis", "pAdj"] <- NA
# 
# write_csv(means, "results/PRJNA285788/immune_gene_identification/gene_count_means.csv")
# 
# gene_counts_imm_ogs <- gene_counts %>% filter(Orthogroup %in% immune_orthogroups)
# 
# osmia_equal <- gene_counts_imm_ogs %>% dplyr::filter(Osmia_bicornis == Osmia_lignaria)
# 
# osmia_higher <- gene_counts_imm_ogs %>% dplyr::filter(Apis_mellifera < Osmia_bicornis &
#   Apis_mellifera < Osmia_lignaria &
#   Bombus_terrestris < Osmia_bicornis &
#   Bombus_terrestris < Osmia_lignaria)
# 
# osmia_lower <- gene_counts_imm_ogs %>% dplyr::filter(Apis_mellifera > Osmia_bicornis &
#   Apis_mellifera > Osmia_lignaria &
#   Bombus_terrestris > Osmia_bicornis &
#   Bombus_terrestris > Osmia_lignaria)
# 
# osmia_higher_vs_drosophila <- osmia_higher %>% 
#   dplyr::filter(Osmia_bicornis > Drosophila_melanogaster &
#                   Osmia_lignaria > Drosophila_melanogaster)
# 
# osmia_lower_vs_drosophila <- osmia_higher %>% 
#   dplyr::filter(Osmia_bicornis < Drosophila_melanogaster &
#                   Osmia_lignaria < Drosophila_melanogaster)
# 
# means_imm_ogs <- as_tibble(apply(gene_counts_imm_ogs[,2:22], MARGIN = 2, mean, na.rm = TRUE)) %>% 
#   mutate(Species = colnames(gene_counts_imm_ogs[,2:22])) %>% dplyr::rename(MeanGeneCount = value) %>% 
#   dplyr::select(Species, MeanGeneCount) %>% 
#   mutate(pVal = apply(gene_counts_imm_ogs[,2:22], MARGIN = 2, wilcox.test_ob)) %>% 
#   mutate(pAdj = p.adjust(pVal)) %>% 
#   arrange(desc(MeanGeneCount))
# 
# means_imm_ogs[means_imm_ogs$Species == "Osmia_bicornis", "pVal"] <- NA
# means_imm_ogs[means_imm_ogs$Species == "Osmia_bicornis", "pAdj"] <- NA
# 
# write_csv(means_imm_ogs, "results/PRJNA285788/immune_gene_identification/gene_count_means_imm_genes.csv")
# 


gene_counts_imm_ogs <- read_tsv("results/PRJNA285788/immune_gene_identification/orthofinder_insects/Phylogenetic_Hierarchical_Orthogroups/N0.tsv") %>% 
  dplyr::filter(OG %in% immune_orthogroups) %>% dplyr::select(-HOG, -"Gene Tree Parent Clade") %>% 
  dplyr::rename(Orthogroup = OG) %>% 
  pivot_longer(cols = Acyrthosiphon_pisum:Vespa_mandarinia, names_to = "Species") %>% 
  dplyr::rename(Proteins = value) %>% mutate(Species = str_replace(Species, "_", " ")) %>% 
  dplyr::filter(Species %in% c("Osmia bicornis", "Osmia lignaria", "Bombus terrestris", 
                               "Apis mellifera", "Drosophila melanogaster"))

gene_counts_imm_ogs_v2 <- gene_counts_imm_ogs %>% mutate(GeneCount = count_names(Proteins)) %>% 
  dplyr::mutate(Category = cut(GeneCount, c(1,2,3,4,1000), labels = c("2", "3", "4", ">4"))) %>% 
  group_by(Species, Category) %>% dplyr::summarise(GeneCount = length(GeneCount)) %>% drop_na() %>% 
  group_by(Species) %>% mutate(GeneCountRel = GeneCount / sum(GeneCount)) 
gene_counts_imm_ogs_v2_text_labels <- gene_counts_imm_ogs_v2 %>% group_by(Species) %>% 
  summarise(GeneCountSum = sum(GeneCount))

gene_counts_imm_ogs_v2.1 <- gene_counts_imm_ogs %>% mutate(GeneCount = count_names(Proteins)) %>% 
  mutate(Category = cut(GeneCount, c(0,1,2,3,1000), labels = c("1", "2", "3", ">3"))) %>% 
  group_by(Species, Category) %>% summarise(GeneCount = length(GeneCount)) %>% drop_na() %>% 
  group_by(Species) %>% dplyr::mutate(GeneCountRel = GeneCount / sum(GeneCount))
gene_counts_imm_ogs_v2.1_text_labels <- gene_counts_imm_ogs_v2.1 %>% group_by(Species) %>% 
  summarise(GeneCountSum = sum(GeneCount))

gene_counts_imm_ogs_v3 <- gene_counts_imm_ogs %>% 
  full_join(dm_gene_group2og) %>% group_by(Species, Group) %>% 
  summarise(GeneCount = join_and_count_unique_names(Proteins)) %>% group_by(Species) %>% 
  summarise(Species, Group, GeneCount,
            AllEqual = ifelse(length(unique(GeneCount)) == 1, TRUE, FALSE)) %>% 
  dplyr::filter(!AllEqual) %>% drop_na() %>% 
  pivot_wider(names_from = Species, values_from = GeneCount) %>% dplyr::select(-AllEqual)
  
rows <- gene_counts_imm_ogs_v3$Group
gene_counts_imm_ogs_v3 <- gene_counts_imm_ogs_v3 %>% dplyr::select(-Group) %>% as.matrix()
row.names(gene_counts_imm_ogs_v3) <- rows

gene_counts_imm_ogs_v3 <- gene_counts_imm_ogs_v3[complete.cases(gene_counts_imm_ogs_v3), ]
gene_counts_imm_ogs_v3_reordered <- gene_counts_imm_ogs_v3[
  hclust(dist(gene_counts_imm_ogs_v3))$order, ] %>% melt() %>% 
  as_tibble() %>% dplyr::rename(Group = Var1, Species = Var2, GeneCount = value)

group_levels <- unique(gene_counts_imm_ogs_v3_reordered$Group)

gene_counts_imm_ogs_no_pw <- gene_counts_imm_ogs_v3_reordered %>% group_by(Group) %>% 
  dplyr::filter(!vectorized_grepl("Pathway", Group), any(GeneCount > 4), 
                !(Group %in% c("ENZYMES"))) %>% 
  drop_na() %>% mutate(Group = factor(Group, levels = group_levels)) %>% 
  mutate(Group = str_to_title(Group)) %>% 
  mutate(Group = replace(Group, Group == "Ras Superfamily Guanine Nucleotide Exchange Factors", 
                         "Ras Guanine Nucl. Exch. Factors")) %>% 
  mutate(Group = replace(Group, Group == "Intracellular Transport Coat, Scaffold And Adaptor Proteins",
                         "Intrac. Transp. Coat, Scaff. & Adapt. Prots.")) %>% 
  mutate(Group = replace(Group, Group == "Chromatin Remodeling Complexes (Atp-Dependent)",
                          "Chromatin Remodeling Complexes"))

gene_counts_imm_ogs_pw <- gene_counts_imm_ogs_v3_reordered %>% group_by(Group) %>% 
  dplyr::filter(vectorized_grepl("Pathway", Group), Group != "BMP Signaling Pathway") %>% 
  drop_na() %>% mutate(Group = factor(Group, levels = group_levels)) %>% 
  mutate(Group = str_remove(Group, " Signaling Pathway"))

colour_palette <- brewer.pal(5, "GnBu")

no_pw_plot <- (ggplot(gene_counts_imm_ogs_no_pw, aes(x = Species, y = Group)) + 
    geom_tile(aes(fill = GeneCount), colour = "white") + 
    geom_text(aes(label = GeneCount)) + 
    geom_vline(xintercept = c(0.5, 5.5)) + 
    scale_fill_gradientn(colours = colour_palette, values = c(0, 0.05, 0.1, 0.3, 0.6, 1),
                         name = "Gene count") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 10, hjust = 1),
          axis.text.y = element_text(size = 10),
          aspect.ratio = 1, axis.title = element_blank(), 
          legend.title = element_text(size = 10), 
          plot.title = element_text(hjust = 0.5)))

pw_plot <- (ggplot(gene_counts_imm_ogs_pw, aes(x = Species, y = Group)) + 
    geom_tile(aes(fill = GeneCount), colour = "white") + 
    geom_text(aes(label = GeneCount)) + 
    geom_vline(xintercept = c(0.5, 5.5)) +  
    scale_fill_gradientn(colours = colour_palette, values = c(0, 0.05, 0.1, 0.3, 0.6, 1),
                         name = "Gene count") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 10, hjust = 1),
          axis.text.y = element_text(size = 10),
          aspect.ratio = 1, axis.title = element_blank(), 
          legend.position = "none", 
          plot.title = element_text(hjust = 0.5)))

barchart <- (ggplot(gene_counts_imm_ogs_v2, aes(x = Species)) + 
    geom_col(aes(y = GeneCountRel, fill = fct_rev(Category))) + 
    scale_fill_hue(direction = -1, name = "Copy number") + 
    scale_x_discrete(labels = c("A.mellifera", "B.terrestris", "D.melanogaster", 
                                "O.bicornis", "O.lignaria")) + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1 )) + 
    ylab("Share of multi-copy orthologues")) + 
    geom_text(data = gene_counts_imm_ogs_v2_text_labels,
              aes(y = 1.04, label = paste(GeneCountSum)))

barchart_v2 <- (ggplot(gene_counts_imm_ogs_v2.1, aes(x = Species)) + 
   geom_col(aes(y = GeneCountRel, fill = fct_rev(Category))) + 
   scale_fill_hue(direction = -1, name = "Copy number") + 
   scale_x_discrete(labels = c("A.mellifera", "B.terrestris", "D.melanogaster", 
                               "O.bicornis", "O.lignaria")) + 
   theme(axis.text.x = element_text(angle = 45, hjust = 1 )) + 
   ylab("Share of orthologues")) +
   geom_text(data = gene_counts_imm_ogs_v2.1_text_labels, 
             aes(y = 1.04, label = paste(GeneCountSum)))


save(pw_plot, no_pw_plot, barchart, 
     file = "results/PRJNA285788/immune_gene_identification/plots/figure2_plots.RData")

ggsave(plot = barchart, "results/PRJNA285788/immune_gene_identification/plots/barchart_imm_genes_no_single_copy_orthl.png",
       width = 25, height = 15, units = "cm")
ggsave(plot = barchart_v2, "results/PRJNA285788/immune_gene_identification/plots/barchart_imm_genes.png",
       width = 25, height = 15, units = "cm")

# diffs <- tibble(Orthgroup = gene_counts$Orthogroup)
# ref_species <- colnames(gene_counts)[2:length(gene_counts)]
# ref_species <- ref_species[!(ref_species %in% c("Osmia_bicornis", "Total"))]
# for(species in ref_species){
#   diffs[[species]] <- gene_counts[[species]] - gene_counts[["Osmia_bicornis"]]
# }
# 
# diffs_long <- diffs %>% pivot_longer(cols = ref_species, names_to = "Species")
# 
# ggplot(diffs_long %>% filter(Species %in% c("Bombus_terrestris", "Apis_mellifera", 
#                                             "Drosophila_melanogaster", "Osmia_lignaria")), 
#        aes(x = value, fill = Species)) + geom_histogram(alpha = 0.5, na.rm = TRUE,
#                                                         position = position_dodge()) + 
#   xlim(-4.1,4.1)
# 
# events <- c("A3SS", "A5SS", "MXE", "RI", "SE")
#   
# splice_variants <- tibble()
# for(event in events){
#   df <- read_tsv(paste0("results/PRJNA285788/sex_differences/diff_splicing/rMATS_results/fromGTF.", 
#                         event, ".txt")) %>% mutate(event = event)
#   splice_variants <- rbind(splice_variants, df %>% dplyr::select(ID, GeneID, event))
# }
# 
# splice_variant_counts <- splice_variants %>% group_by(GeneID) %>% summarise(VariantCount = length(GeneID)) %>% 
#   arrange(desc(VariantCount))
# 
# splice_variant_counts_imm_genes <- splice_variant_counts %>% filter(GeneID %in% immune_genes)
# splice_variant_counts_non_imm_genes <- splice_variant_counts %>% filter(!(GeneID %in% immune_genes))
# 
# print(mean(splice_variant_counts_imm_genes$VariantCount))
# print(mean(splice_variant_counts_non_imm_genes$VariantCount))
# print(wilcox.test(splice_variant_counts_imm_genes$VariantCount,
#             splice_variant_counts_non_imm_genes$VariantCount))
# 
# splice_variant_counts_event_lvl <- splice_variants %>% group_by(GeneID, event) %>% summarise(VariantCount = length(GeneID)) %>% 
#   arrange(desc(VariantCount))
# 
# splice_variant_counts_event_lvl_imm_genes <- splice_variant_counts_event_lvl %>% 
#   filter(GeneID %in% immune_genes)
# splice_variant_counts_event_lvl_non_imm_genes <- splice_variant_counts_event_lvl %>% 
#   filter(!(GeneID %in% immune_genes))
# 
# for(Event in events){
#   a <- splice_variant_counts_event_lvl_imm_genes %>% filter(event == Event) %>% pull(VariantCount)
#   b <- splice_variant_counts_event_lvl_non_imm_genes %>% filter(event == Event) %>% pull(VariantCount)
#   print(Event)
#   print(mean(a))
#   print(mean(b))
#   print(wilcox.test(a,b))
# }

gene_counts <- read_tsv("results/PRJNA285788/immune_gene_identification/orthofinder_insects/Phylogenetic_Hierarchical_Orthogroups/N0.tsv") %>% 
  dplyr::select(-HOG, -"Gene Tree Parent Clade") %>% dplyr::rename(Orthogroup = OG) %>% 
  pivot_longer(cols = Acyrthosiphon_pisum:Vespa_mandarinia, names_to = "Species") %>% 
  dplyr::rename(Proteins = value) %>% mutate(Species = str_replace(Species, "_", " ")) %>% 
  dplyr::filter(Species %in% c("Osmia bicornis", "Osmia lignaria", "Bombus terrestris", 
                               "Apis mellifera", "Drosophila melanogaster"))

single_copy_orthologues <-read_tsv("results/PRJNA285788/immune_gene_identification/orthofinder_bees/Phylogenetic_Hierarchical_Orthogroups/N0.tsv") %>% 
  dplyr::select(-HOG, -"Gene Tree Parent Clade") %>% dplyr::rename(Orthogroup = OG) %>% 
  pivot_longer(cols = Apis_mellifera:Osmia_lignaria, names_to = "Species") %>% 
  dplyr::rename(Proteins = value) %>% mutate(Species = str_replace(Species, "_", " ")) %>% 
  group_by(Species, Orthogroup) %>% summarise(GeneCount = join_and_count_unique_names(Proteins)) %>% 
  group_by(Orthogroup) %>% dplyr::filter(all(GeneCount == 1)) %>% pull(Orthogroup) %>% unique()

write_lines(x = single_copy_orthologues, 
            file = "results/PRJNA285788/immune_gene_identification/immune_genes/single_copy_orthologues.txt")

gene_counts_v2 <- gene_counts %>% mutate(GeneCount = count_names(Proteins)) %>% 
  dplyr::mutate(Category = cut(GeneCount, c(1,2,3,4,1000), labels = c("2", "3", "4", ">4"))) %>% 
  group_by(Species, Category) %>% summarise(GeneCount = length(GeneCount)) %>% drop_na() %>% 
  group_by(Species) %>% dplyr::mutate(GeneCountRel = GeneCount / sum(GeneCount)) 

gene_counts_v3 <- gene_counts %>% 
  full_join(dm_gene_group2og) %>% group_by(Species, Group) %>% 
  summarise(GeneCount = join_and_count_unique_names(Proteins)) %>% group_by(Species) %>% 
  summarise(Species, Group, GeneCount,
            AllEqual = ifelse(length(unique(GeneCount)) == 1, TRUE, FALSE)) %>% 
  dplyr::filter(!AllEqual) %>% drop_na() %>% 
  pivot_wider(names_from = Species, values_from = GeneCount) %>% dplyr::select(-AllEqual)

rows <- gene_counts_v3$Group
gene_counts_v3 <- gene_counts_v3 %>% dplyr::select(-Group) %>% as.matrix()
row.names(gene_counts_v3) <- rows

gene_counts_v3 <- gene_counts_v3[complete.cases(gene_counts_v3), ]
gene_counts_v3_reordered <- gene_counts_v3[
  hclust(dist(gene_counts_v3))$order, ] %>% melt() %>% 
  as_tibble() %>% dplyr::rename(Group = Var1, Species = Var2, GeneCount = value)

group_levels <- unique(gene_counts_v3_reordered$Group)

gene_counts_no_pw <- gene_counts_v3_reordered %>% group_by(Group) %>% 
  dplyr::filter(!vectorized_grepl("Pathway", Group), any(GeneCount > 4), 
                !(Group %in% c("ENZYMES"))) %>% 
  drop_na() %>% mutate(Group = factor(Group, levels = group_levels)) %>% 
  mutate(Group = str_to_title(Group)) %>% 
  mutate(Group = replace(Group, Group == "Ras Superfamily Guanine Nucleotide Exchange Factors", 
                         "Ras Guanine Nucl. Exch. Factors")) %>% 
  mutate(Group = replace(Group, Group == "Intracellular Transport Coat, Scaffold And Adaptor Proteins",
                         "Intrac. Transp. Coat, Scaff. & Adapt. Prots.")) %>% 
  mutate(Group = replace(Group, Group == "Chromatin Remodeling Complexes (Atp-Dependent)",
                         "Chromatin Remodeling Complexes"))

gene_counts_pw <- gene_counts_v3_reordered %>% group_by(Group) %>% 
  dplyr::filter(vectorized_grepl("Pathway", Group), Group != "BMP Signaling Pathway") %>% 
  drop_na() %>% mutate(Group = factor(Group, levels = group_levels)) %>% 
  mutate(Group = str_remove(Group, " Signaling Pathway"))

colour_palette <- brewer.pal(5, "GnBu")

no_pw_plot_all_genes <- (ggplot(gene_counts_no_pw, aes(x = Species, y = Group)) + 
                 geom_tile(aes(fill = GeneCount), colour = "white") + 
                 geom_text(aes(label = GeneCount)) + 
                 geom_vline(xintercept = c(0.5, 5.5)) + 
                 scale_fill_gradientn(colours = colour_palette, values = c(0, 0.05, 0.1, 0.3, 0.6, 1),
                                      name = "Gene count") +
                 theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 10, hjust = 1),
                       axis.text.y = element_text(size = 10),
                       aspect.ratio = 1, axis.title = element_blank(), 
                       legend.title = element_text(size = 10), 
                       plot.title = element_text(hjust = 0.5)))

pw_plot_all_genes <- (ggplot(gene_counts_pw, aes(x = Species, y = Group)) + 
              geom_tile(aes(fill = GeneCount), colour = "white") + 
              geom_text(aes(label = GeneCount)) + 
              geom_vline(xintercept = c(0.5, 5.5)) +  
              scale_fill_gradientn(colours = colour_palette, values = c(0, 0.05, 0.1, 0.3, 0.6, 1),
                                   name = "Gene count") +
              theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 10, hjust = 1),
                    axis.text.y = element_text(size = 10),
                    aspect.ratio = 1, axis.title = element_blank(), 
                    legend.position = "none", 
                    plot.title = element_text(hjust = 0.5)))

barchart_all_genes <- (ggplot(gene_counts_v2, aes(x = Species, y = GeneCountRel, fill = fct_rev(Category))) + 
               geom_col() + 
               scale_fill_hue(direction = -1, name = "Copy number") + 
               scale_x_discrete(labels = c("A.mellifera", "B.terrestris", "D.melanogaster", 
                                           "O.bicornis", "O.lignaria")) + 
               theme(axis.text.x = element_text(angle = 45, hjust = 1 )) + 
               ylab("Share of multi-copy orthogroups"))

save(pw_plot, no_pw_plot, barchart, 
     file = "results/PRJNA285788/immune_gene_identification/plots/figure2_plots_all_genes.RData")


gene_counts_imm_ogs_v4 <- gene_counts_imm_ogs %>% 
  full_join(dm_gene_group2og) %>% 
  group_by(Orthogroup, Species) %>% summarise(Orthogroup, Species, Group, 
                                     GeneCounts = ifelse(any(!is.na(Proteins)), 1, 0)) %>% 
  group_by(Species, Group) %>% summarise(GeneCounts = sum(GeneCounts)) %>% group_by(Species) %>% 
  summarise(Species, Group, GeneCounts,
            AllEqual = ifelse(length(unique(GeneCounts)) == 1, TRUE, FALSE)) %>% 
  dplyr::filter(!AllEqual) %>% drop_na() %>% 
  mutate(GeneCounts = as.numeric(str_remove(GeneCounts, " "))) %>% 
  pivot_wider(names_from = Species, values_from = GeneCounts) %>% dplyr::select(-AllEqual) %>% 
  ungroup()


rows <- gene_counts_imm_ogs_v4$Group
gene_counts_imm_ogs_v4 <- gene_counts_imm_ogs_v4 %>% dplyr::select(-Group) %>% as.matrix()
row.names(gene_counts_imm_ogs_v4) <- rows

gene_counts_imm_ogs_v4 <- gene_counts_imm_ogs_v4[complete.cases(gene_counts_imm_ogs_v4), ]
gene_counts_imm_ogs_v4_reordered <- gene_counts_imm_ogs_v4[
  hclust(dist(gene_counts_imm_ogs_v4))$order, ] %>% melt() %>% 
  as_tibble() %>% dplyr::rename(Group = Var1, Species = Var2, GeneCount = value)

group_levels <- unique(gene_counts_imm_ogs_v4_reordered$Group)

gene_counts_imm_ogs_no_pw_v2 <- gene_counts_imm_ogs_v4_reordered %>% group_by(Group) %>% 
  dplyr::filter(!vectorized_grepl("Pathway", Group), any(GeneCount > 4), 
                !(Group %in% c("ENZYMES"))) %>% 
  drop_na() %>% mutate(Group = factor(Group, levels = group_levels)) %>% 
  mutate(Group = str_to_title(Group)) %>% 
  mutate(Group = replace(Group, Group == "Ras Superfamily Guanine Nucleotide Exchange Factors", 
                         "Ras Guanine Nucl. Exch. Factors")) %>% 
  mutate(Group = replace(Group, Group == "Intracellular Transport Coat, Scaffold And Adaptor Proteins",
                         "Intrac. Transp. Coat, Scaff. & Adapt. Prots.")) %>% 
  mutate(Group = replace(Group, Group == "Chromatin Remodeling Complexes (Atp-Dependent)",
                         "Chromatin Remodeling Complexes"))

gene_counts_imm_ogs_pw_v2 <- gene_counts_imm_ogs_v4_reordered %>% group_by(Group) %>% 
  dplyr::filter(vectorized_grepl("Pathway", Group), Group != "BMP Signaling Pathway") %>% 
  drop_na() %>% mutate(Group = factor(Group, levels = group_levels)) %>% 
  mutate(Group = str_remove(Group, " Signaling Pathway"))

colour_palette <- brewer.pal(5, "GnBu")

no_pw_plot_v2 <- (ggplot(gene_counts_imm_ogs_no_pw_v2, aes(x = Species, y = Group)) + 
                 geom_tile(aes(fill = GeneCount), colour = "white") + 
                 geom_text(aes(label = GeneCount)) + 
                 geom_vline(xintercept = c(0.5, 5.5)) + 
                 scale_fill_gradientn(colours = colour_palette, values = c(0, 0.05, 0.1, 0.3, 0.6, 1),
                                      name = "Gene count") +
                 theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 10, hjust = 1),
                       axis.text.y = element_text(size = 10),
                       aspect.ratio = 1, axis.title = element_blank(), 
                       legend.title = element_text(size = 10), 
                       plot.title = element_text(hjust = 0.5)))

pw_plot_v2 <- (ggplot(gene_counts_imm_ogs_pw_v2, aes(x = Species, y = Group)) + 
              geom_tile(aes(fill = GeneCount), colour = "white") + 
              geom_text(aes(label = GeneCount)) + 
              geom_vline(xintercept = c(0.5, 5.5)) +  
              scale_fill_gradientn(colours = colour_palette, values = c(0, 0.05, 0.1, 0.3, 0.6, 1),
                                   name = "Gene count") +
              theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 10, hjust = 1),
                    axis.text.y = element_text(size = 10),
                    aspect.ratio = 1, axis.title = element_blank(), 
                    legend.position = "none", 
                    plot.title = element_text(hjust = 0.5)))
