#' Find all manifested descendants of a node
#' @export
#' @param node node in \code{tree}
#' @param tree \code{phylo} class object
#' @param nodebank output from \code{\link{nodeBank}}
#' @examples 
#' set.seed(1)
#' library(ape)
#' tree <- rtree(10)
#' nb <- nodeBank(11,tree,0.8)
#' manifestChildren(11,tree,nb)
manifestChildren <- function(node,tree,nodebank){
  m <- nodebank[node==node & manifest==T,node]
  if (length(m)==0){
    ch <- phangorn::Descendants(tree,node,'children')
    while(length(ch)>0){
      m <- c(m,nodebank[node %in% ch & manifest==T,node])
      ch <- nodebank[node %in% ch & manifest==F,node]
      if (length(ch)>0){
        ch <- lapply(ch,FUN=function(nd,tr) phangorn::Descendants(tree,nd,'children')) %>%
          unlist
        ch <- setdiff(ch,1:length(tree$tip.label))
      }
    }
  }
  return(m)
}