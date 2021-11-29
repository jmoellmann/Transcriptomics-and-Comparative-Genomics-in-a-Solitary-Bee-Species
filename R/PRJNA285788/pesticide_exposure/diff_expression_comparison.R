library(dplyr)
library(DESeq2)
library(eulerr)
library(ggplot2)

project_dir <- "/media/jannik/Samsung_T5/masters_project/ApisTranscriptomics/"
setwd(project_dir)

dds_kallisto <- readRDS("results/PRJNA285788/pesticide_exposure/diff_expr/kallisto/dds_wald.RDS")
dds_star <- readRDS("results/PRJNA285788/pesticide_exposure/diff_expr/STAR/dds_wald.RDS")

immune_genes <- read.table(
  "results/PRJNA285788/immune_gene_identification/immune_genes/ob_merged_imm_genes.txt")[["V1"]]

# PCAs ----------------------------------------------------------------------------------------

pca_kallisto <- plotPCA(vst(dds_kallisto), intgroup = "Group")
pca_star <- plotPCA(vst(dds_star), intgroup = "Group")

ggsave("results/PRJNA285788/pesticide_exposure/diff_expr/kallisto/pca.png", plot = pca_kallisto)
ggsave("results/PRJNA285788/pesticide_exposure/diff_expr/STAR/pca.png", plot = pca_star)

# Pesticide treatment stats -------------------------------------------------------------
heatmaps <- list()

for(treatment in c("Imidacloprid", "Thiacloprid")){
  res_kallisto <- results(dds_kallisto, contrast = c("Group", "Control", treatment))
  res_star <- results(dds_star, contrast = c("Group", "Control", treatment))
  
  res_genes <- res_star %>% as_tibble() %>% mutate(GeneID = row.names(res_star)) %>% 
    dplyr::select(GeneID, baseMean:padj)
  
  write.csv(as.data.frame(res_kallisto), paste0(
    "results/PRJNA285788/pesticide_exposure/diff_expr/kallisto/res_", treatment, ".csv"))
  write.csv(as.data.frame(res_star), paste0(
    "results/PRJNA285788/pesticide_exposure/diff_expr/STAR/res_", treatment, ".csv"))

  # DEG overlap
  degs_overlap <- list()
  degs_overlap[["kallisto"]] <- row.names(res_kallisto[!is.na(res_kallisto$padj) & 
                                                         res_kallisto$padj < 0.05, ])
  degs_overlap[["star"]] <- row.names(res_star[!is.na(res_star$padj) & res_star$padj < 0.05, ])
  
  deg_overlap_plot <- plot(eulerr::euler(degs_overlap), quantities = TRUE, 
                           main = "Overlap in DEGs (padj < 0.05)")
  
  # Immune DEG overlap
  immune_degs_overlap <- list()
  immune_degs_overlap[["kallisto"]] <- row.names(res_kallisto[!is.na(res_kallisto$padj) & 
                                                                res_kallisto$padj < 0.05 & 
                                                      row.names(res_kallisto) %in% immune_genes, ])
  immune_degs_overlap[["star"]] <- row.names(res_star[!is.na(res_star$padj) & 
                                                        res_star$padj < 0.05 &
                                                        row.names(res_star) %in% immune_genes, ])
  
  immune_deg_overlap_plot <- plot(eulerr::euler(immune_degs_overlap), quantities = TRUE, 
                                  main = "Overlap in immune DEGs (padj < 0.05)")
  
  ggsave(plot = deg_overlap_plot, filename = 
           paste0("results/PRJNA285788/pesticide_exposure/diff_expr/deg_overlap_", treatment, ".png"), 
         device = "png")
  ggsave(plot = immune_deg_overlap_plot, filename = 
           paste0("results/PRJNA285788/pesticide_exposure/diff_expr/imm_deg_overlap_", 
                  treatment, ".png"), 
         device = "png")
  
  
  # Immune DEG stats
  
  res_degs <- res_genes %>% dplyr::filter(padj < 0.05)
  
  res_immune_genes <- res_genes %>% dplyr::filter(GeneID %in% immune_genes)
  res_immune_degs <- res_genes %>% dplyr::filter(GeneID %in% immune_genes, padj < 0.05)

  
  # Heatmap -------------------------------------------------------------------------------------
  
  if(treatment == "Imidacloprid"){
    samples <- c(paste0("SRR", 7286081:7286088)) 
  } else {
    samples <- c(paste0("SRR", c(7286081:7286084, 7286089:7286092))) 
  }
  
  
  vst_counts <- DESeq2::vst(DESeq2::counts(dds_star))
  
  counts <- vst_counts[row.names(vst_counts) %in% immune_genes, samples]
  
  normalised_counts <- t(scale(t(counts), center = TRUE, scale = TRUE))
  normalised_counts <- normalised_counts[complete.cases(normalised_counts), ]
  normalised_counts_reordered <- normalised_counts[hclust(dist(normalised_counts))$order, ]
  
  normalised_counts_tbl <- melt(normalised_counts_reordered) %>% 
    as_tibble() %>% 
    dplyr::rename("GeneID" = "Var1", "Sample" = "Var2", "ReadCount" = "value") %>% 
    dplyr::filter(GeneID %in% res_immune_degs$GeneID)
  
  heatmaps[[treatment]] <- (ggplot(normalised_counts_tbl, aes(Sample, GeneID)) + 
                              geom_tile(aes(fill = ReadCount), colour = "white") + 
                              scale_fill_gradient(name = "Normalised \nRead Counts") + 
                              scale_x_discrete(labels = c(paste0("Control ", 1:4), paste0(treatment, 1:4))) + 
                              theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 10, hjust = 1),
                                    axis.text.y = element_blank(),
                                    aspect.ratio = 1, axis.title = element_blank(),
                                    axis.ticks.y = element_blank()))
  
  ggsave(paste0("results/PRJNA285788/pesticide_exposure/diff_expr/deg_heatmap_", 
                tolower(treatment), ".png"), heatmaps[[treatment]],
         width = 25, height = 15, unit = "cm")
  
  
  
  test_matrix <- matrix(data = c(length(res_immune_genes$GeneID), length(res_immune_degs$GeneID), 
                                 length(res_genes$GeneID), length(res_degs$GeneID)), 
                        nrow = 2, ncol = 2)
  f.test <- fisher.test(test_matrix)
  
  out_file <- paste0("results/PRJNA285788/pesticide_exposure/diff_expr/stats_", treatment, ".txt")
  sink(out_file, append = FALSE)
  cat(paste("Genes:", test_matrix[1,2], "\n"))
  cat(paste("Immune Genes:", test_matrix[1,1], "\n"))
  cat(paste("DEGs:", test_matrix[2,2], "\n"))
  cat(paste("Immune DEGs:", test_matrix[2,1], "\n"))
  cat("-------------------------------------------------", "\n")
  cat(paste("DEGs up in control group:", 
            nrow(res_degs[res_degs$log2FoldChange > 0, ]), "\n"))
  cat(paste("DEGs up in pesticide treated group:", 
            nrow(res_degs[res_degs$log2FoldChange < 0, ]), "\n"))
  cat(paste("Immune DEGs up in control group:", 
            nrow(res_immune_degs[res_immune_degs$log2FoldChange > 0, ]), "\n"))
  cat(paste("Immune DEGs up in pesticide treated group:", 
            nrow(res_immune_degs[res_immune_degs$log2FoldChange < 0, ]), "\n"))
  cat("-------------------------------------------------", "\n")
  cat(paste0("Share of immune genes diff. expr.: ", 
             round((test_matrix[2,1] / test_matrix[1,1]), 4) * 100, "%", "\n"))
  cat(paste0("Share of all genes diff. expr.: ", 
             round((test_matrix[2,2] / test_matrix[1,2]), 4) * 100, "%", "\n"))
  cat("-------------------------------------------------", "\n")
  cat("\n")
  cat(paste0("Fisher-test p-value: ", round(f.test$p.value, 3), "\n"))
  sink()
  
}

save(heatmaps, file = "results/PRJNA285788/pesticide_exposure/diff_expr/deg_heatmaps.RData")
