# Libraries -----------------------------------------------------------------------------------
library(biomaRt)
library(tidyr)
library(stringr)
library(readr)
library(eulerr)
library(ggplot2)
library(dplyr)

# Functions -----------------------------------------------------------------------------------
gather_coords <- function(df){
  coords <- c()
  for(i in 1:nrow(df)){
    if(df[i, "event"] == "SE"){
      tmp <- paste0(df[i, "upstreamEE"], "-", df[i, "exonStart_0base"] + 1, ":",
                    df[i, "exonEnd"], "-", df[i, "downstreamES"] + 1)
    } else if(df[i, "event"] == "RI"){
      tmp <- paste0(df[i, "upstreamES"] + 1, "-", df[i, "upstreamEE"], ":", 
                    df[i, "downstreamES"] + 1, "-", df[i, "downstreamEE"])
    } else if(df[i, "event"] == "MXE"){
      tmp <- paste0(df[i, "upstreamEE"], "-", df[i, "X1stExonStart_0base"] + 1, ":",
                    df[i, "X1stExonEnd"], "-", df[i, "downstreamES"] + 1, ":", 
                    df[i, "upstreamEE"], "-", df[i, "X2ndExonStart_0base"] + 1, ":", 
                    df[i, "X2ndExonEnd"], "-", df[i, "downstreamES"] + 1)
    } else if(df[i, "event"] == "A3SS" && df[i, "strand"] == "-" || 
              df[i, "event"] == "A5SS" && df[i, "strand"] == "+"){
      tmp <-paste0(df[i, "longExonEnd"], "-", df[i, "flankingES"] + 1, ":", 
                   df[i, "shortEE"], "-", df[i, "flankingES"] + 1)
    } else if(df[i, "event"] == "A3SS" && df[i, "strand"] == "+" ||
              df[i, "event"] == "A5SS" && df[i, "strand"] == "-"){
      tmp <- paste0(df[i, "flankingEE"], "-", df[i, "longExonStart_0base"] + 1, ":", 
                    df[i, "flankingEE"], "-", df[i, "shortES"] + 1)
    }
    coords <- c(coords, tmp)
  }
  return(coords)
}

getSplicingResults <- function(data_dir){
  results <- list()
  body_parts <- c("Thorax", "Abdomen", "Head")
  groups <- c("adult_male", "virgin_queen", "mated_queen", "worker_forager", "worker_nurse")
  splicing_events <- c("SE", "RI", "MXE", "A3SS", "A5SS")
  for(body_part in body_parts){
    for(i in 1:(length(groups)-1)){
      for(j in (i+1):length(groups)){
        for(splicing_event in splicing_events){
          groups_name <- paste0(groups[i], "_", groups[j])
          res_table <- as_tibble(read.table(paste(data_dir, "rMATS", body_part, groups_name, 
                                                  paste0(splicing_event, ".MATS.JCEC.txt"), 
                                                  sep = "/"), header = TRUE)) %>% 
            mutate(event = splicing_event) %>%
            mutate(coords = gather_coords(.)) %>% 
            rename("p.val" = "PValue", "dPSI" = "IncLevelDifference", 
                   "cond1" = "IncLevel1", "cond2" = "IncLevel2") %>% 
            mutate(event = replace(event, event == "A3SS", "A3")) %>% 
            mutate(event = replace(event, event == "A5SS", "A5")) %>% 
            mutate(event = replace(event, event == "MXE", "MX")) %>% 
            dplyr::select(GeneID, chr, coords, strand, event, p.val, FDR, cond1, cond2, dPSI) %>% 
            arrange(coords)
          res_table <- res_table %>% mutate(GeneID = str_remove(GeneID, "gene:"))
          results[[body_part]][[groups_name]][[splicing_event]] <- res_table
        }
      }
    }
  }
  return(results)
}

# Set paths -----------------------------------------------------------------------------------
project_dir <- "/media/jannik/Samsung_T5/masters_project/ApisTranscriptomics"
kallisto_dir <- paste0(project_dir, "/results/PRJNA533756/after_trimming/diff_splicing/kallisto")
star_dir <- paste0(project_dir, "/results/PRJNA533756/after_trimming/diff_splicing/STAR")

# Start analysis ------------------------------------------------------------------------------------
results <- list()

# Kallisto ------------------------------------------------------------------------------------
results[["kallisto"]] <- getSplicingResults(data_dir = kallisto_dir)

# STAR ----------------------------------------------------------------------------------------
results[["STAR"]] <- getSplicingResults(data_dir = star_dir)

# Save results --------------------------------------------------------------------------------
saveRDS(results[["kallisto"]], paste0(project_dir, "/results/PRJNA533756/after_trimming/diff_splicing/kallisto/rMATS/res.RDS"))
saveRDS(results[["STAR"]], paste0(project_dir, "/results/PRJNA533756/after_trimming/diff_splicing/STAR/rMATS/res.RDS"))
