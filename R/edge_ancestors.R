#' get path of ancestral edges all the way to, but not including, root
#' @export
#' @param edge integer - which edge
#' @param tree \code{phylo} class object
#' @examples
#' 
edge_ancestors <- function(edge,tree){
  rt=ape::Ntip(tree)+1
  anc=tree$edge[edge,1]
  if (anc==rt){
    return(NULL)
  } else {
    anc_edges <- NULL
    while(!anc==rt){
      new_anc_edge <- which(tree$edge[,2]==anc)
      anc_edges <- c(anc_edges,new_anc_edge)
      anc=tree$edge[new_anc_edge,1]
    }
    return(anc_edges)
  }
}