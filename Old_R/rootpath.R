#' find path from node to root
#' @export
#' @param node node of \code{tree}
#' @param tree phylo class object
#' @examples 
#' library(dendromap)
#' set.seed(1)
#' tree <- rtree(20)
#' rp=rootpath(33,tree)
#' 
#' plot(tree)
#' nodelabels(rp,rp)
rootpath <- function(node,tree){
  N=length(tree$tip.label)
  nds <- node
  i=1
  while(!(N+1) %in% nds){
    ix=tree$edge[,2]==nds[i]
    if (!any(ix)){
      stop(paste('Could not find find descendant node',nds[i],'in tree$edge'))
    } else {
      nds <- c(nds,tree$edge[ix,1])
      i=i+1
    }
  }
  return(nds)
}
