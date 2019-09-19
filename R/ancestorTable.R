#' get table of ancestral node-pair scores
#' @export
#' @param U score matrix
#' @param row.node node in \code{row.tree}
#' @param col.node node in \code{col.tree}
#' @param row.tree phylo class object whose i'th node (starting from the root being 1) corresponds to i'th row of \code{U}
#' @param col.tree phylo class object whose j'th node (starting from the root being 1) corresponds to j'th column of \code{U}
ancestorTable <- function(U, row.node, col.node, row.tree, col.tree){
  m <- length(row.tree$tip.label)
  n <- length(col.tree$tip.label)
  orientation <- function(parent_node,daughter_node,tree){
    ix <- which(tree$edge[,1]==parent_node)
    if (daughter_node==tree$edge[ix[1],2]){
      return(1)
    } else if (daughter_node==tree$edge[ix[2],2]){
      return(-1)
    } else {
      return(0)
    }
  }
  if (row.node==(m+1) | col.node==(n+1)){
    return(NULL)
  } else {
    row.nodepath <- as.character(nodepath(row.tree, from = row.node, to=m+1))
    row.nodepath <- as.numeric(row.nodepath[2:(length(row.nodepath))])
    
    col.nodepath <- as.character(nodepath(col.tree, from = col.node, to=n+1))
    col.nodepath <- as.numeric(col.nodepath[2:(length(col.nodepath))])

    ancestor.table <- U[row.nodepath-m,col.nodepath-n]
    
    n.row <- length(row.nodepath)
    n.col <- length(col.nodepath)
    
    row.orientations <- mapply(FUN=orientation,row.nodepath,
                               c(row.node,row.nodepath[1:(n.row-1)]),
                               MoreArgs = list('tree'=row.tree))
    names(row.orientations) <- row.nodepath
    col.orientations <- mapply(FUN=orientation,col.nodepath,
                               c(col.node,col.nodepath[1:(n.col-1)]),
                               MoreArgs = list('tree'=col.tree))
    names(col.orientations) <- col.nodepath
    orientation.matrix <- t(col.orientations %*% t(row.orientations))
    colnames(orientation.matrix) <- colnames(ancestor.table)
    rownames(orientation.matrix) <- rownames(ancestor.table)
    return(list('table'=ancestor.table,'orientation'=orientation.matrix))
  }
}
