#' Get scores for ancestor sequences
#' @export
#' @param U score matrix whose i'th row corresponds to the (starting from root=1) i'th node in \code{row.tree} and j'th column corresponds to the j'th node in \code{col.tree}
#' @param row.node node in \code{row.tree}
#' @param col.node node in \code{col.tree}
#' @param row.tree phylo class object
#' @param col.tree phylo class object
#' @param orientation.match Logical. Whether orientation of paths must be perfect match. Currently only supports \code{TRUE} until soft-matching is coded up.
#' @examples
#' library(dendromap)
#' set.seed(1)
#' row.tree <- ape::rtree(50)
#' col.tree <- ape::rtree(10)
#' U <- t(treeBasis(row.tree)) %*% matrix(rnorm(500),nrow=50) %*% treeBasis(col.tree)
#' par(mfrow=c(1,2))
#' plot(row.tree)
#' nodelabels()
#' nodelabels(77,77,bg='red')
#' plot(col.tree)
#' nodelabels()
#' nodelabels(16,16,bg='red')
#' 
#' Anc <- matchAncestors(U,77,16,row.tree,col.tree)

matchAncestors <- function(U,row.node,col.node,row.tree,col.tree,orientation.match=TRUE){
  
  Anc <- ancestorTable(U,row.node,col.node,row.tree,col.tree)
  ancestor_table <- Anc$table
  ancestor_orientations <- Anc$orientation
  rm('Anc')
  if (orientation.match){
    ancestor_table[ancestor_table*ancestor_orientations<0] <- 0
    rm('ancestor_orientations')
    gc()
  }
  
  if (all(ancestor_table==0)){
    return(NULL)
  } else {
    
    rownames(ancestor_table) <- gsub('node_','',rownames(ancestor_table))
    colnames(ancestor_table) <- gsub('node_','',colnames(ancestor_table))
    best_scoring_words <- list()
    row_node = c()
    col_node = c()
    score = c()
    if (is.null(row.tree$bootstrap.values)){
      row.tree$bootstrap.values <- rep(1,ape::Nnode(row.tree))
    }
    
    max_word_length <- min(ncol(ancestor_table),nrow(ancestor_table))
    
    for (word_length in 1:max_word_length){
      combinations <- combn(colnames(ancestor_table), word_length)
      
      for (comb in 1:ncol(combinations)){
        comb_name <- combinations[,comb]
        comb_score <- word_length
        
        index_max = 1
        row_matches = c()
        col_matches = c()
        
        for (i in comb_name){
          
          row_match <- rownames(ancestor_table)[index_max:nrow(ancestor_table)][which.max(ancestor_table[,as.character(i)])]
          col_match <- as.character(i)
          row_matches = c(row_matches, row_match)
          col_matches = c(col_matches, col_match)
          scores <- ancestor_table[,col_match]
          
          valid_scores <- intersect(which(scores!=0),index_max:length(scores))
          
          if (length(valid_scores)>0){
            best_ix <- which.max(abs(scores[valid_scores]))
            best_score <- scores[valid_scores][best_ix]
            nd <- as.numeric(rownames(ancestor_table)[valid_scores[best_ix]])
            node_bootstrap <- row.tree$bootstrap_values[nd-length(row.tree$tip.label)]
            
            ## Can change line below to prod for quantiles
            # comb_score <- prod(c(comb_score, best_score, (node_bootstrap/100)))
            comb_score <- comb_score+best_score
        
            ##
            
            index_max = which.max(ancestor_table[,col_match][valid_scores]) -1
          } else {
            comb_score <- 0
          }
        }
        row_node = c(row_node, paste(row_matches, collapse = '--'))
        col_node = c(col_node, paste(col_matches, collapse = '--'))
        score = c(score, comb_score)
      }
    }
    data = data.table(row_node, col_node, score)
    data$N <- lapply(data$row_node,strsplit,'--') %>% lapply(unlist) %>% sapply(length)
    return(data)
  }
}
