library(tidyverse)
library(ggplot2)
library(ggpubr)

project_dir <- "/run/media/T5/masters_project/01_main_analysis/"
setwd(project_dir)

get_merged_df <- function(in_path, in_file_name){
   merged_df <- tibble()
   for (ontology in c("CC", "MF", "BP")){
      if(ontology == "BP"){
         n <- 20
      } else {
         n <- 5
      }
      merged_df <- read_tsv(paste0(in_path, ontology, "/", in_file_name)) %>% 
         arrange(weight01KS) %>% mutate(ontology = ontology) %>% dplyr::slice(1:n) %>% 
         rbind(merged_df)
   }
   return(merged_df)   
}

stretch_terms <- function(terms, lim = 50){
   stretched_terms <- c()
   for (term in terms) {
      n <- nchar(term)
      if(n < lim){
         stretched_term <- paste0(str_flatten(rep(" ", 50 - n)), term)
         stretched_terms <- c(stretched_terms, stretched_term)
      } else {
         stretched_terms <- c(stretched_terms, term)
      }   
   }
   return(stretched_terms)
}


# Differential expression ---------------------------------------------------------------------

male_up_df <- get_merged_df("04_results/PRJNA285788/sex_differences/diff_expr/STAR/funct_enrich/", 
                            "gsea_dm_ns20_lfc.tsv") %>%
   mutate(Term = factor(Term, levels = rev(Term)))
female_up_df <- get_merged_df("04_results/PRJNA285788/sex_differences/diff_expr/STAR/funct_enrich/", 
                              "gsea_dm_ns20_lfc_neg.tsv") %>%
   mutate(Term = ifelse(nchar(Term) > 50, paste0(str_sub(Term, 1, 50), "[..]"), Term)) %>% 
   mutate(Term = factor(Term, levels = rev(Term)))

male_up_go_plot <- ggplot(male_up_df, aes(x = Term, y = -log2(weight01KS), fill = ontology)) + 
   geom_col() + 
   xlab("GO term") + ylab("-log2(p)") + 
   geom_hline(yintercept = c(-log(0.05, 2)), linetype = "dashed") + 
   coord_flip() #+ 
#   scale_y_discrete(limits = rev(levels(TermAsFactor)))

female_up_go_plot <- ggplot(female_up_df, aes(x = Term, y = -log2(weight01KS), fill = ontology)) + 
   geom_col() + 
   xlab("GO term") + ylab("-log2(p)") + 
   geom_hline(yintercept = c(-log(0.05, 2)), linetype = "dashed") + 
   scale_x_discrete(expand = c(0, 0)) +    
   coord_flip() #+ 
#   scale_y_discrete(limits = rev(levels(TermAsFactor)))

ggsave("04_results/PRJNA285788/sex_differences/diff_expr/STAR/funct_enrich/male_up.png", 
       male_up_go_plot)

ggsave("04_results/PRJNA285788/sex_differences/diff_expr/STAR/funct_enrich/female_up.png", 
       female_up_go_plot)

combined_plot <- ggpubr::ggarrange(male_up_go_plot, female_up_go_plot, nrow = 2, labels = "AUTO")

ggsave("04_results/PRJNA285788/sex_differences/diff_expr/STAR/funct_enrich/male_female_combined.png",
       combined_plot, width = 8.22, height = 7.0)


# Differential splicing -----------------------------------------------------------------------

splicing_df <- get_merged_df("04_results/PRJNA285788/sex_differences/diff_splicing/funct_enrich/",
              "gsea_dm_ns20.tsv") %>% 
   mutate(Term = factor(Term, levels = rev(Term)))

splicing_go_plot <- ggplot(splicing_df, aes(x = Term, y = -log2(weight01KS), fill = ontology)) + 
   geom_col() + 
   xlab("GO term") + ylab("-log2(p)") + 
   geom_hline(yintercept = c(-log(0.05, 2)), linetype = "dashed") + 
   coord_flip() #+ 


ggsave("04_results/PRJNA285788/sex_differences/diff_splicing/funct_enrich/splicing.png", 
       splicing_go_plot)


# Pesticide exposure differential expression --------------------------------------------------


# Differential expression ---------------------------------------------------------------------

for (pesticide in c("thiacloprid", "imidacloprid")) {
   pesticide_down_df <- get_merged_df(paste0(
      "04_results/PRJNA285788/pesticide_exposure/diff_expr/STAR/funct_enrich/", pesticide, "/"), 
      "gsea_dm_ns20_lfc_pos.tsv") %>% mutate(Term = factor(Term, levels = rev(Term)))
   pesticide_up_df <- get_merged_df(paste0(
      "04_results/PRJNA285788/pesticide_exposure/diff_expr/STAR/funct_enrich/", pesticide, "/"),
      "gsea_dm_ns20_lfc_neg.tsv") %>%
      #mutate(Term = ifelse(nchar(Term) > 50, paste0(str_sub(Term, 1, 50), "[..]"), Term)) %>% 
      mutate(Term = factor(Term, levels = rev(Term)))
   
   pesticide_down_go_plot <- ggplot(pesticide_down_df, aes(x = Term, y = -log2(weight01KS), 
                                                       fill = ontology)) + 
      geom_col() + 
      xlab("GO term") + ylab("-log2(p)") + 
      geom_hline(yintercept = c(-log(0.05, 2)), linetype = "dashed") + 
      coord_flip() #+ 
   #   scale_y_discrete(limits = rev(levels(TermAsFactor)))
   
   pesticide_up_go_plot <- ggplot(pesticide_up_df, aes(x = Term, y = -log2(weight01KS), 
                                                           fill = ontology)) + 
      geom_col() + 
      xlab("GO term") + ylab("-log2(p)") + 
      geom_hline(yintercept = c(-log(0.05, 2)), linetype = "dashed") + 
      scale_x_discrete(expand = c(0, 0)) +    
      coord_flip() #+ 
   #   scale_y_discrete(limits = rev(levels(TermAsFactor)))
   
   ggsave(paste0(
      "04_results/PRJNA285788/pesticide_exposure/diff_expr/STAR/funct_enrich/", pesticide, 
      "/pesticide_down.png"), pesticide_down_go_plot, width = 8.22, height = 5)
   
   ggsave(paste0(
      "04_results/PRJNA285788/pesticide_exposure/diff_expr/STAR/funct_enrich/", pesticide, 
      "/pesticide_up.png"), pesticide_up_go_plot, width = 8.22, height = 5)
   
   pesticide_combined_plot <- ggpubr::ggarrange(pesticide_up_go_plot, pesticide_down_go_plot,
                                                nrow = 2, labels = "AUTO")
   
   ggsave(paste0(
      "04_results/PRJNA285788/pesticide_exposure/diff_expr/STAR/funct_enrich/", pesticide, 
          "/combined.png"), pesticide_combined_plot, width = 8.22, height = 7.0)
}


# Differential splicing -----------------------------------------------------------------------

for (pesticide in c("thiacloprid", "imidacloprid")) {
   pesticide_splicing_df <- get_merged_df(paste0("04_results/PRJNA285788/pesticide_exposure/diff_splicing/", pesticide, 
                                                 "/funct_enrich/"), "gsea_dm_ns20.tsv") %>% 
      mutate(Term = factor(Term, levels = rev(Term)))
   
   pesticide_splicing_go_plot <- ggplot(pesticide_splicing_df, aes(x = Term, y = -log2(weight01KS), fill = ontology)) + 
      geom_col() + 
      xlab("GO term") + ylab("-log2(p)") + 
      geom_hline(yintercept = c(-log(0.05, 2)), linetype = "dashed") + 
      coord_flip() #+ 
   
   
   ggsave(paste0("04_results/PRJNA285788/pesticide_exposure/diff_splicing/", pesticide, 
                 "/funct_enrich/splicing.png"), pesticide_splicing_go_plot)
}



# Immune gene analysis ------------------------------------------------------------------------


