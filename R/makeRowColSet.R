#' make RowCol list for parallelized input to findNodeSeq
#' @param U matrix of bilinear mappings from {row.tree}x{col.tree} nodes
#' @param forbidden.nodes integer vector of row nodes to exclude from consideration
makeRowColSet <- function(U,forbidden.nodes=NULL){
  nds <- as.numeric(sapply(rownames(U),FUN=function(x) strsplit(x,'_')[[1]][2]))
  col.nds <- as.numeric(sapply(colnames(U),FUN=function(x) strsplit(x,'_')[[1]][2]))
  forbidden.nodes <- c(forbidden.nodes,nds[rowSums(U==0)==ncol(U)])
  rows <- setdiff(1:nrow(U),which(nds %in% forbidden.nodes))
  cols <- setdiff(1:ncol(U),col.nds[colSums(U==0)==nrow(U)])
  RCset <- expand.grid(rows,cols)
    split(x = RCset,seq(nrow(RCset))) %>% lapply(unlist) %>%
    return()
}
