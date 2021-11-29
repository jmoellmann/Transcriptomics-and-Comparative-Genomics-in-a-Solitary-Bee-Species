library(topGO)

project_dir <- "/run/media/T5/masters_project/01_main_analysis/"

source(paste0(project_dir, "02_scripts/general/GO_mappings_to_list.R"))

get_gsea_results <- function(stats, gene2GO, mode = "ks", alpha = 0.05, padj = FALSE,
                             nodeSize = 20, ontology = "BP", outputTopGOobject = FALSE){
  
  if(mode == "fisher"){
    topGOobject <- new("topGOdata", ontology = ontology,
                       allGenes = stats,
                       nodeSize = nodeSize, annot = annFUN.gene2GO, gene2GO = gene2GO)
    
    resultsFisher <- runTest(topGOobject, algorithm = "classic", statistic = "fisher")
    #resultsFisher.weight01 <- runTest(topGOobject, algorithm = "weight01", statistic = "fisher")
    
    topNodes = length(score(resultsFisher))
    
    allRes <- GenTable(topGOobject, 
                       classicFisher = resultsFisher, orderBy = "classicFisher", 
                       ranksOf = "classicFisher", topNodes = topNodes, numChar = 500)
    
    allRes$classicFisher <- as.numeric(str_remove(allRes$classicFisher, "< "))
    #allRes$weight01FisherPadj = p.adjust(allRes$weigh01Fisher, method = "BH")
    if(padj) {
       allRes$classicFisherPadj <- p.adjust(allRes$classicFisher, method = "BH")
    }
    
  } else if(mode == "ks"){
    if(alpha == 0){
      topGOobject <- new("topGOdata", ontology = ontology,
                         allGenes = stats, geneSel = function(x) x,
                         nodeSize = nodeSize, annot = annFUN.gene2GO, gene2GO = gene2GO)
      } else {
      topGOobject <- new("topGOdata", ontology = ontology,
                         allGenes = stats, geneSel = function(x) x < alpha,
                         nodeSize = nodeSize, annot = annFUN.gene2GO, gene2GO = gene2GO)
      }
      #resultsKS <- runTest(topGOobject, algorithm = "classic", statistic = "ks")
      resultsKS.weight01 <- runTest(topGOobject, algorithm = "weight01", statistic = "ks")
      
      topNodes = length(score(resultsKS.weight01))
      
      allRes <- GenTable(topGOobject, weight01KS = resultsKS.weight01, orderBy = "weight01KS", 
                         ranksOf = "weight01KS", topNodes = topNodes, numChar = 500)
      
      allRes$weight01KS <- as.numeric(str_remove(allRes$weight01KS, "< "))
      if(padj){
         allRes$weight01KSpadj = p.adjust(allRes$weight01KS, method = "BH")
         #allRes$classicKSpadj = p.adjust(allRes$classicKS, method = "BH")
      }
      
  } else if(mode == "both"){
     if(alpha == 0){
        topGOobject <- new("topGOdata", ontology = ontology,
                           allGenes = stats, geneSel = function(x) x,
                           nodeSize = nodeSize, annot = annFUN.gene2GO, gene2GO = gene2GO)
     } else {
        topGOobject <- new("topGOdata", ontology = ontology,
                           allGenes = stats, geneSel = function(x) x < alpha,
                           nodeSize = nodeSize, annot = annFUN.gene2GO, gene2GO = gene2GO)
     }
     
     #resultsKS <- runTest(topGOobject, algorithm = "classic", statistic = "ks")
     resultsKS.weight01 <- runTest(topGOobject, algorithm = "weight01", statistic = "ks")
     resultsFisher <- runTest(topGOobject, algorithm = "classic", statistic = "fisher")
     
     topNodes = length(score(resultsKS.weight01))
     
     allRes <- GenTable(topGOobject, weight01KS = resultsKS.weight01, classicFisher = resultsFisher,
                        orderBy = "weight01KS", ranksOf = "weight01KS", topNodes = topNodes, 
                        numChar = 500)
     
     allRes$classicFisher <- as.numeric(str_remove(allRes$classicFisher, "< "))
     allRes$weight01KS <- as.numeric(str_remove(allRes$weight01KS, "< "))
     if(padj){
        allRes$weight01KSpadj = p.adjust(allRes$weight01KS, method = "fdr")
        allRes$classicFisherpadj = p.adjust(allRes$classicFisher, method = "fdr")
        #allRes$classicKSpadj = p.adjust(allRes$classicKS, method = "BH")
     }
  }
  if(outputTopGOobject){
    return(list("results" = allRes, "object" = topGOobject))
  } else {
    return(allRes)
    }
}
