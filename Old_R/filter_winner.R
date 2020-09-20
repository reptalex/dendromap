#' filter winner from lineages found
#' @export
#' @param winner index - element of Lineages maximizing score
#' @param Lineages list of lineages to filter
#' @param row.tree phylo object
#' @param rc_table \code{\link{makeRCtable}}
#' @param row.nodemap \code{\link{makeNodeMap}}
filter_winner <- function(winner,Lineages,row.tree,
                          rc_table,row.nodemap){
  row.node <- rc_table[rc_index %in% Lineages[[winner]],min(row.node)]
  rt <- min(row.nodemap$node)
  Desc <- row.node:(row.nodemap[node==row.node,row.node+pos+neg])
  anc <- ape::nodepath(row.tree,row.node,rt)
  
  filtered_rc_ix <- rc_table[row.node %in% unique(c(Desc,anc)),rc_index]
  
  Lineages <- sapply(Lineages,FUN=function(a,b) setdiff(a,b),b=filtered_rc_ix)
  Lineages <- Lineages[sapply(Lineages,FUN=function(x) length(x)>0)]
  return(Lineages)
}
