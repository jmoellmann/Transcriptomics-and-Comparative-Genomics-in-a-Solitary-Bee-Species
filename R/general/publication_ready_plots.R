library(ggplot2)
library(cowplot)
library(ggpubr)
library(tidyverse)

setwd("/media/jannik/Samsung_T5/masters_project/ApisTranscriptomics/")


# Figure 1 ------------------------------------------------------------------------------------

load("results/PRJNA285788/immune_gene_identification/plots/figure1_plots.RData")

(figure_1 <- cowplot::plot_grid(overlap_imm_genes_plot_no_title, dm_stats2ob_stats_plot_v2, dm_stats2ob_stats_plot_v1, 
           unique_imm_gene_dm_seq_ident_distr_plot, scale = c(0.9,0.95,0.95,0.95), 
           labels = c("A", "B", "C", "D"), rel_widths = c(1,1.09)))

# figure_1_v2 <- ggarrange(overlap_imm_genes_plot_no_title, dm_stats2ob_stats_plot_v2, dm_stats2ob_stats_plot_v1, 
#                          unique_imm_gene_dm_seq_ident_distr_plot)

ggsave(plot = figure_1, filename = "results/PRJNA285788/immune_gene_identification/plots/figure1.png",
       width = 25, height = 15, units = "cm", device = "png", bg = "white")


# Figure 2 ------------------------------------------------------------------------------------

load("results/PRJNA285788/immune_gene_identification/plots/figure2_plots.RData")


(figure_2 <- cowplot::plot_grid(pw_plot, no_pw_plot, rel_widths = c(0.68, 1), 
                                label_y = c(1, 1), label_x = c(0.05, 0.05), 
                                scale = c(0.95, 0.95), labels = c("A", "B")))

ggsave(plot = figure_2, filename = "results/PRJNA285788/immune_gene_identification/plots/figure2.png",
       width = 34, height = 15, units = "cm", device = "png", bg = "white")


# Figure 3 ------------------------------------------------------------------------------------

load("results/PRJNA285788/sex_differences/diff_expr/heatmap.RData")

sex_diff_go_df <- read_tsv("results/PRJNA285788/sex_differences/diff_expr/STAR/funct_enrich/BP/gsea_imm_gene_ns20.tsv") %>% 
  mutate(stat = -log(weight01KS, 2)) %>% mutate(Term = paste0(Term, " (", Annotated, ")")) %>% 
  dplyr::filter(weight01KS < 0.045) %>% arrange(stat) %>% 
  mutate(TermAsFactor = factor(Term, levels = Term))

sex_diff_go_plot <- ggplot(sex_diff_go_df, aes(x = TermAsFactor, y = stat)) + geom_col(fill = "steelblue1", col = "steelblue4") + xlab("GO Term") +
  ylab("-log2(p)") + geom_hline(yintercept = c(-log(0.05, 2)), linetype = "dashed") + coord_flip()

(figure_3 <- cowplot::plot_grid(heatmap, sex_diff_go_plot, rel_widths = c(0.8, 1),
                                label_y = c(1, 1), label_x = c(0.05, 0.05), 
                                scale = c(0.90, 0.90), labels = c("A", "B")))

ggsave(plot = figure_3, filename = "results/PRJNA285788/sex_differences/diff_expr/figure3.png",
       width = 35, height = 15, units = "cm", device = "png", bg = "white")

# Figure 4 ------------------------------------------------------------------------------------

load("results/PRJNA285788/pesticide_exposure/diff_expr/deg_heatmaps.RData")

imid_exp_go_df <- read_tsv("results/PRJNA285788/pesticide_exposure/diff_expr/STAR/funct_enrich/imidacloprid/BP/gsea_imm_gene_ns50.tsv") %>% 
  mutate(stat = -log(weight01KS, 2)) %>% mutate(Term = paste0(Term, " (", Annotated, ")")) %>% 
  dplyr::filter(weight01KS < 0.045) %>% arrange(stat) %>% 
  mutate(TermAsFactor = factor(Term, levels = Term))

imid_exp_go_plot <- ggplot(imid_exp_go_df, aes(x = TermAsFactor, y = stat)) + geom_col(fill = "steelblue1", col = "steelblue4") + xlab("GO Term") +
  ylab("-log2(p)") + geom_hline(yintercept = c(-log(0.05, 2)), linetype = "dashed") + coord_flip()

thia_exp_go_df <- read_tsv("results/PRJNA285788/pesticide_exposure/diff_expr/STAR/funct_enrich/thiacloprid/BP/gsea_imm_gene_ns50.tsv") %>% 
  mutate(stat = -log(weight01KS, 2)) %>% mutate(Term = paste0(Term, " (", Annotated, ")")) %>% 
  dplyr::filter(weight01KS < 0.045) %>% arrange(stat) %>% 
  mutate(TermAsFactor = factor(Term, levels = Term))

thia_exp_go_plot <- ggplot(thia_exp_go_df, aes(x = TermAsFactor, y = stat)) + geom_col(fill = "steelblue1", col = "steelblue4") + xlab("GO Term") +
  ylab("-log2(p)") + geom_hline(yintercept = c(-log(0.05, 2)), linetype = "dashed") + coord_flip()

(figure_4 <- cowplot::plot_grid(heatmaps[["Imidacloprid"]], imid_exp_go_plot, 
                                heatmaps[["Thiacloprid"]], thia_exp_go_plot, 
                                rel_widths = c(0.7, 1),
                                label_y = c(1, 1), label_x = c(0.05, 0.05), 
                                scale = c(0.90, 0.90, 0.90, 0.90), labels = c("A", "B", "C", "D")))

ggsave(plot = figure_4, filename = "results/PRJNA285788/pesticide_exposure/diff_expr/figure4.png",
       width = 30, height = 18, units = "cm", device = "png", bg = "white")
