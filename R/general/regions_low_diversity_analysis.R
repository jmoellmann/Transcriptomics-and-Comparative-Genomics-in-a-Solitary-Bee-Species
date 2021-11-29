#library(pafr)
library(ggplot2)
#library(ggpubr)
library(tidyverse)
library(gggenes)

setwd("/run/media/T5/masters_project/01_main_analysis/")
# 
# ob_ol <- pafr::read_paf("01_data/02_phylogenetic_analysis/chr_rearrangements/ob_ol_genome_alignment.paf")
# 
# ob_ol_long_ali <- subset(ob_ol, alen > 1e4 & mapq > 40)
# prim_ob_ol_ali <- filter_secondary_alignments(ob_ol_long_ali)
# 
# ggplot(prim_ob_ol_ali, aes(alen, dv)) + 
#    geom_point(alpha=0.6, colour="steelblue", size=2) + 
#    scale_x_continuous("Alignment length (kb)", label =  function(x) x/ 1e3) +
#    scale_y_continuous("Per base divergence") + 
#    theme_pubr()
# 
# by_q <- aggregate(dv ~ qname, data=prim_ob_ol_ali, FUN=mean)
# knitr::kable(by_q)
# 
# 
# bt_am <- pafr::read_paf("01_data/02_phylogenetic_analysis/chr_rearrangements/bt_am_genome_alignment.paf")
# 
# bt_am_long_ali <- subset(bt_am, alen > 1e4 & mapq > 40)
# prim_bt_am_ali <- filter_secondary_alignments(bt_am_long_ali)
# 
# ggplot(prim_bt_am_ali, aes(alen, dv)) + 
#    geom_point(alpha=0.6, colour="steelblue", size=2) + 
#    scale_x_continuous("Alignment length (kb)", label =  function(x) x/ 1e3) +
#    scale_y_continuous("Per base divergence") + 
#    theme_pubr()
# 
# by_q <- aggregate(dv ~ qname, data=prim_bt_am_ali, FUN=mean)
# knitr::kable(by_q)
# 
# bt_am_rold <- pafr::read_paf("01_data/02_phylogenetic_analysis/chr_rearrangements/bt_am_rold_alignment.paf")
# 
# bt_am_rold_long_ali <- subset(bt_am_rold, alen > 1e4 & mapq > 40)
# prim_bt_am_rold_ali <- filter_secondary_alignments(bt_am_rold_long_ali)
# 
# ggplot(prim_bt_am_rold_ali, aes(alen, dv)) + 
#    geom_point(alpha=0.6, colour="steelblue", size=2) + 
#    scale_x_continuous("Alignment length (kb)", label =  function(x) x/ 1e3) +
#    scale_y_continuous("Per base divergence") + 
#    theme_pubr()
# 
# by_q <- aggregate(dv ~ qname, data=prim_bt_am_ali, FUN=mean)
# knitr::kable(by_q)
# 
# 
# bt_ob_rold <- pafr::read_paf("01_data/02_phylogenetic_analysis/chr_rearrangements/bt_ob_rold_alignment.paf")
# 
# bt_ob_rold_long_ali <- subset(bt_ob_rold, alen > 1e4 & mapq > 40)
# prim_bt_ob_rold_ali <- filter_secondary_alignments(bt_ob_rold_long_ali)
# 
# ggplot(prim_bt_ob_rold_ali, aes(alen, dv)) + 
#    geom_point(alpha=0.6, colour="steelblue", size=2) + 
#    scale_x_continuous("Alignment length (kb)", label =  function(x) x/ 1e3) +
#    scale_y_continuous("Per base divergence") + 
#    theme_pubr()
# 
# by_q <- aggregate(dv ~ qname, data=prim_bt_ob_ali, FUN=mean)
# knitr::kable(by_q)
# 
# plot_synteny(bt_ob_rold, q_chrom = "OU015516.1", t_chrom = "NC_015762.1:350000-570000", centre = TRUE) + theme_bw()
# 

# Low diversity region in Osmia ---------------------------------------------------------------

genes <- readRDS("04_results/PRJNA285788/genes.RDS")
single_copy_ogs <- read_table("01_data/03_reference/Orthogroups_SingleCopyOrthologues.txt")[[1]]


merged_rold_df <- bind_rows(genes$Bombus_terrestris %>% filter(lowDivGene, OG != "NA") %>% 
                               mutate(species = "B. terrestris"),
                            genes$Osmia_bicornis %>% filter(lowDivGene, OG != "NA") %>% 
                               mutate(species = "O. bicornis")) %>% 
   mutate(species = factor(species, levels = c("B. terrestris", "O. bicornis")))

merged_rold_df %>% ggplot(aes(xmin = start, xmax = end, y = species, fill = OG, 
                              forward = ifelse(strand == "+", TRUE, FALSE))) + 
   geom_gene_arrow(arrowhead_height = unit(3, "mm"), arrowhead_width = unit(0, "mm"))+ 
   xlim(c(360000, 520000)) + theme_genes() + theme(legend.position = "none")
