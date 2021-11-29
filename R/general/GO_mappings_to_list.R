# Function To transform R data frame with duplicated entries in the first column to R list ----

# requires the GO ID mapping to be input as data frame without empty cells
# with the gene IDs in the first column and the GO IDs in the second column
GO_mappings_to_list <- function(df){
  results_list <-  list()
  j = 1
  k = 1
  go_terms <- c(df[1, 2])
  gene_names <- c()
  for(i in 2:nrow(df)){
    if(df[i, 1] == df[j, 1]){
      go_terms <- c(go_terms, df[i, 2])
    } else {
      results_list[[k]] <- go_terms
      gene_names <- c(gene_names, df[j, 1])
      k = k + 1
      j = i
      go_terms <- c(df[i, 2])
    }
  }
  results_list[[k]] <- go_terms
  gene_names <- c(gene_names, df[j, 1])
  names(results_list) <- gene_names
  
  return(results_list)
}