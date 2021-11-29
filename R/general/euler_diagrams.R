library(eulerr)
library(DESeq2)

setwd("/home/jannik/mogon/")

# Euler plot with diff. expressed genes -------------------------------------------------------

adult_deseq <- results(readRDS("results/GSE85433_drone_worker/after_trimming/diff_expr/STAR/deseq2_ds.rds"))
pupa_deseq <- results(readRDS("results/GSE45408_pupa/after_trimming/diff_expr/STAR/deseq2_ds.rds"))
larva_deseq <- results(readRDS("results/GSE61253_larva/after_trimming/diff_expr/STAR/deseq2_ds.rds"))

adult_deseq_up <-row.names(adult_deseq[(adult_deseq$log2FoldChange > 0) & (!is.na(adult_deseq$pvalue)) & 
                                (adult_deseq$pvalue < 0.05), ])
pupa_deseq_up <-row.names(pupa_deseq[(pupa_deseq$log2FoldChange > 0) & (!is.na(pupa_deseq$pvalue)) & 
                                         (pupa_deseq$pvalue < 0.05), ])
larva_deseq_up <-row.names(larva_deseq[(larva_deseq$log2FoldChange > 0) & (!is.na(larva_deseq$pvalue)) & 
                                         (larva_deseq$pvalue < 0.05), ])

adult_deseq_down <-row.names(adult_deseq[(adult_deseq$log2FoldChange < 0) & (!is.na(adult_deseq$pvalue)) & 
                                           (adult_deseq$pvalue < 0.05), ])
pupa_deseq_down <-row.names(pupa_deseq[(pupa_deseq$log2FoldChange < 0) & (!is.na(pupa_deseq$pvalue)) & 
                                           (pupa_deseq$pvalue < 0.05), ])
larva_deseq_down <-row.names(larva_deseq[(larva_deseq$log2FoldChange < 0) & (!is.na(larva_deseq$pvalue)) & 
                                           (larva_deseq$pvalue < 0.05), ])

adult_deseq <- row.names(adult_deseq[!is.na(adult_deseq$pvalue) & adult_deseq$pvalue < 0.05, ])
pupa_deseq <- row.names(pupa_deseq[!is.na(pupa_deseq$pvalue) & pupa_deseq$pvalue < 0.05, ])
larva_deseq <- row.names(larva_deseq[!is.na(larva_deseq$pvalue) & larva_deseq$pvalue < 0.05, ])


euler_deseq <- plot(euler(list(adult = adult_deseq, pupa = pupa_deseq, larva = larva_deseq)), 
                    quantities = TRUE)

euler_deseq_up <- plot(euler(list(adult = adult_deseq_up, pupa = pupa_deseq_up, 
                                  larva = larva_deseq_up)), quantities = TRUE)

euler_deseq_down <- plot(euler(list(adult = adult_deseq_down, pupa = pupa_deseq_down, 
                                    larva = larva_deseq_down)), quantities = TRUE)

# Euler plot with GO terms --------------------------------------------------------------------

adult_goterms <- read.csv("results/GSE85433_drone_worker/after_trimming/diff_expr/STAR/go_enrichment/GO_BP_dm_orthologues.csv")
pupa_goterms <- read.csv("results/GSE45408_pupa/after_trimming/diff_expr/STAR/go_enrichment/GO_BP_dm_orthologues.csv")
larva_goterms <- read.csv("results/GSE61253_larva/after_trimming/diff_expr/STAR/go_enrichment/GO_BP_dm_orthologues.csv")

euler_goterms <- plot(euler(list(adult = adult_goterms$GO.ID, pupa = pupa_goterms$GO.ID, 
                                 larva = larva_goterms$GO.ID)), quantities = TRUE)

dir.create("results/general/diff_expr/STAR/euler", showWarnings = FALSE, recursive = TRUE)

png("results/general/diff_expr/STAR/euler/euler_deseq.png", width = 800, height = 600)
euler_deseq
dev.off()

png("results/general/diff_expr/STAR/euler/euler_deseq_up.png", width = 800, height = 600)
euler_deseq_up
dev.off()

png("results/general/diff_expr/STAR/euler/euler_deseq_down.png", width = 800, height = 600)
euler_deseq_down
dev.off()

png("results/general/diff_expr/STAR/euler/euler_goterms.png", width = 800, height = 600)
euler_goterms
dev.off()
