library(tidyverse)
library(gggenes)
library(biomaRt)
library(topGO)


apply_ks_test <- function(x, y){
   if(length(na.omit(x)) > 0) { 
      return(ks.test(x, y)[["p.value"]]) 
   } else {
      return(1)
   }
}

get_sliding_window_stats <- function(genes, chr_len, window_size, window_shift){
   
   read_counts <- c()
   chrs <- c()
   starts <- c()
   n_genes_tot <- c()
   n_degs_tot <- c()
   mean_read_counts <- c()
   median_read_counts <- c()
   pvals_tot <- list()
   raw_read_counts <- list()
   for (curr_chr in chr_len$chr[str_starts(chr_len$chr, "OU0")]) {
      i = 0
      max_len = chr_len %>% filter(chr == curr_chr) %>% pull(len)
      print(paste0(curr_chr, ": ", max_len))
      while ((i + window_size) < max_len){
         curr_genes <- genes %>% filter(chr == curr_chr, start >= i, start < i + window_size)
         n_genes <- length(curr_genes %>% pull(meanReadCount))
         n_degs <- curr_genes %>% filter(padj < 0.05) %>% nrow()
         median_read_count <- median(curr_genes %>% pull(meanReadCount), na.rm = TRUE)
         mean_read_count <- mean(curr_genes %>% pull(meanReadCount), na.rm = TRUE)
         read_counts_raw <- curr_genes %>% pull(meanReadCount)
         pvals <- curr_genes %>% pull(padj)
         starts <- c(starts, i)
         chrs <- c(chrs, curr_chr)
         n_genes_tot <- c(n_genes_tot, n_genes)
         n_degs_tot <- c(n_degs_tot, n_degs)
         mean_read_counts <- c(mean_read_counts, mean_read_count)
         median_read_counts <- c(median_read_counts, median_read_count)
         raw_read_counts[[paste0(curr_chr, "-", as.character(i))]] <- read_counts_raw
         pvals_tot[[paste0(curr_chr, "-", as.character(i))]] <- pvals
         i = i + window_shift
      }
   }
   res <- tibble(read_counts = raw_read_counts, chr = chrs, start = starts, 
                 n_genes = n_genes_tot, n_degs = n_degs_tot, 
                 mean_read_counts = mean_read_counts, 
                 median_read_counts = median_read_counts, pvals = pvals_tot)

   read_count_pvals <- sapply(res$read_counts, function(x) apply_ks_test(x, genes$meanReadCount))
   read_count_padj_vals <- p.adjust(read_count_pvals, method = "fdr")
   deg_pvals <- sapply(res$pvals, function(x) apply_ks_test(x, genes$pvalue))
   deg_padj_vals <- p.adjust(deg_pvals)
   res$read_count_pvals <- read_count_pvals
   res$deg_pvals <- deg_pvals
   res$read_count_padj_vals <- read_count_padj_vals
   res$deg_padj_vals <- deg_padj_vals

   return(res)
}


get_groups <- function(column){
   i = 1
   groups <- c(1)
   diffs <- diff(column)
   for(diff in diffs){
      if(abs(diff) > 200000){
         i = i +1
      }
      groups <- c(groups, i)
   }
   return(groups)
}


apply_wilcox_test <- function(test_data){
   out_vector <- c()
   for (i in 1:length(test_data)) {
      if(i %% 500 == 0){ print(paste0("Applying Mann-Whitney-U-test. Processed: ", i)) }
      out_vector <- c(out_vector, wilcox.test(mean(test_data), test_data[i])[["p.value"]])
   }
   return(out_vector)
}

project_dir <- "/run/media/T5/masters_project/01_main_analysis/"
hyphy_res_dir <- paste0(project_dir, "01_data/02_phylogenetic_analysis/hyphy_results/")

setwd(project_dir)
source("02_scripts/PRJNA285788/general/get_gsea_results.R")

genes <- readRDS("04_results/PRJNA285788/genes.RDS")
single_copy_ogs <- read_table("01_data/03_reference/Orthogroups_SingleCopyOrthologues.txt")[[1]]

# Comparison of species -----------------------------------------------------------------------

osmia_shared_selected_ogs <- intersect(genes$Osmia_lignaria %>% filter(episDivSelection) %>% pull(OG), 
                                       genes$Osmia_bicornis %>% filter(episDivSelection) %>% pull(OG))

# Regions of low/high expression (sex data) ----------------------------------------------------


ob_chr_len <- read_tsv("01_data/03_reference/annotation/ob_chr_lengths.txt", 
                       col_names = c("chr", "len"))

res <- get_sliding_window_stats(genes$Osmia_bicornis, chr_len = ob_chr_len, 
                                window_size = 1000000, window_shift = 200000)

pval_vs_ngenes_plot <- res %>% ggplot(aes(x = n_genes, y = read_count_pvals)) + geom_point() + 
   geom_smooth() + xlab("Gene count") + ylab("P-value") + 
   geom_hline(aes(yintercept = 0.05), linetype = "dashed")

ggsave("04_results/PRJNA285788/pval_vs_ngenes.png", pval_vs_ngenes_plot)

res$chr_v2 <- res$chr
res$chr <- factor(str_replace_all(res$chr, 
                                  c("OU015504.1" = "Chr 1", "OU015505.1" = "Chr 2",
                                    "OU015506.1" = "Chr 3", "OU015507.1" = "Chr 4",
                                    "OU015508.1" = "Chr 5", "OU015509.1" = "Chr 6",
                                    "OU015510.1" = "Chr 7", "OU015511.1" = "Chr 8",
                                    "OU015512.1" = "Chr 9", "OU015513.1" = "Chr 10",
                                    "OU015514.1" = "Chr 11", "OU015515.1" = "Chr 12",
                                    "OU015516.1" = "Chr 13", "OU015517.1" = "Chr 14",
                                    "OU015518.1" = "Chr 15", "OU015519.1" = "Chr 16")),
                  levels = c(paste0("Chr ", 1:16)))

genes$Osmia_bicornis %>% filter(str_detect(chr, "OU0")) %>%
   ggplot(aes(x = start, y = log2(meanReadCount), color = ifelse(padj < 0.05, TRUE, FALSE))) +
   geom_point(alpha = 0.7) + facet_wrap(vars(chr)) +
   scale_color_discrete(name = "DEG")

all_chr_plot <- res %>% filter(str_detect(chr, "Chr")) %>% 
   ggplot(aes(x = start, y = median_read_counts, 
              color = ifelse(read_count_padj_vals < 0.01, TRUE, FALSE))) + 
   geom_point(alpha = 0) + geom_spoke(aes(angle = 0, radius = 1e6), size = 1.5, alpha = 0.7) + 
   scale_x_continuous(breaks = c(0, 5e6, 1e7, 1.5e7), labels = c("0Mb", "5Mb", "10Mb", "15Mb")) + 
   facet_wrap(vars(chr)) + scale_color_discrete(name = "FDR < 0.01\nKS test") + 
   ylab("Average read count") + xlab("Chromosome position") # +

   #ggtitle("Under-/overexpressed regions on all chromosome")

ggsave("04_results/PRJNA285788/all_chr_plot.png", all_chr_plot)

genomic_regions <- res %>% filter(read_count_padj_vals < 0.01) %>% mutate(group = get_groups(start)) %>% 
   group_by(group) %>% summarise(chr, median_exp = median(median_read_counts), 
                                 start = min(start), n_genes = sum(n_genes), n = n()) %>% 
   distinct() %>% mutate(length = 1e6 + (n - 1) * 2e5, end = start + length) %>% ungroup() %>%
   dplyr::select(chr, start, end, length, n_genes, median_exp) %>% 
   mutate(median_exp = round(median_exp, digits = 2)) %>% arrange(chr, start)

med_read_count <- median(genes$Osmia_bicornis$meanReadCount)

genomic_upexpressed_regions <- genomic_regions %>% filter(median_exp > med_read_count)
genomic_downexpressed_regions <- genomic_regions %>% filter(median_exp < med_read_count)

genomic_upexpressed_regions_plot <- gridExtra::tableGrob(
   genomic_upexpressed_regions, cols = c("Chromosome", "Start", "End", "Length", 
                                         "Gene count", "Med. expr."))
genomic_downexpressed_regions_plot <- gridExtra::tableGrob(
   genomic_downexpressed_regions, cols = c("Chromosome", "Start", "End", "Length", 
                                          "Gene count", "Med. expr."))

ggsave("04_results/PRJNA285788/upexpressed_regions.png", genomic_upexpressed_regions_plot, 
       width = 7.23, height = 6.23)
ggsave("04_results/PRJNA285788/downexpressed_regions.png", genomic_downexpressed_regions_plot, 
       width = 7.23, height = 7.73)

write_csv(genomic_upexpressed_regions, "04_results/PRJNA285788/upexpressed_regions.csv")
write_csv(genomic_downexpressed_regions, "04_results/PRJNA285788/downexpressed_regions.csv")


# Regions of low / high expression (pesticide data) -------------------------------------------

ob_genes_thiacloprid <- readRDS("04_results/PRJNA285788/thiacloprid_genes.RDS")
ob_genes_imidacloprid <- readRDS("04_results/PRJNA285788/imidacloprid_genes.RDS")

thiacloprid_res <- get_sliding_window_stats(ob_genes_thiacloprid, chr_len = ob_chr_len, 
                                window_size = 1000000, window_shift = 200000)

thiacloprid_res$chr_v2 <- thiacloprid_res$chr
thiacloprid_res$chr <- factor(str_replace_all(thiacloprid_res$chr, 
                                  c("OU015504.1" = "Chr 1", "OU015505.1" = "Chr 2",
                                    "OU015506.1" = "Chr 3", "OU015507.1" = "Chr 4",
                                    "OU015508.1" = "Chr 5", "OU015509.1" = "Chr 6",
                                    "OU015510.1" = "Chr 7", "OU015511.1" = "Chr 8",
                                    "OU015512.1" = "Chr 9", "OU015513.1" = "Chr 10",
                                    "OU015514.1" = "Chr 11", "OU015515.1" = "Chr 12",
                                    "OU015516.1" = "Chr 13", "OU015517.1" = "Chr 14",
                                    "OU015518.1" = "Chr 15", "OU015519.1" = "Chr 16")),
                  levels = c(paste0("Chr ", 1:16)))

thiacloprid_all_chr_plot <- thiacloprid_res %>% filter(str_detect(chr, "Chr")) %>% 
   ggplot(aes(x = start, y = median_read_counts, 
              color = ifelse(read_count_padj_vals < 0.01, TRUE, FALSE))) + 
   geom_point(alpha = 0) + geom_spoke(aes(angle = 0, radius = 1e6), size = 1.5, alpha = 0.7) + 
   scale_x_continuous(breaks = c(0, 5e6, 1e7, 1.5e7), labels = c("0Mb", "5Mb", "10Mb", "15Mb")) + 
   facet_wrap(vars(chr)) + scale_color_discrete(name = "FDR < 0.01\nKS test") + 
   ylab("Average read count") + xlab("Chromosome position") # +

ggsave("04_results/PRJNA285788/thiacloprid_all_chr_plot.png", thiacloprid_all_chr_plot)
