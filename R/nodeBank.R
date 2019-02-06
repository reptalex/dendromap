#' Make a dataset of nodes in a tree and whether/not they are manifested
#' 
#' @export
#' @param node node in tree
#' @param tree \code{phylo} class object
#' @param prob probability between 0 and 1
#' @output data frame with columns "node", "manifest" and "terminal"
#' @examples
#' set.seed(1)
#' library(ape)
#' tree <- rtree(10)
#' nb <- nodeBank(11,tree,0.5)
nodeBank <- function(node,tree,prob=1){
  nb <- data.table('node'=phangorn::Descendants(tree,node,'all'))
  nb <- nb[node %in% (ape::Ntip(tree)+1:tree$Nnode)]
  nb[,manifest:=sample(c(T,F),size=.N,replace = T,prob=c(prob,1-prob))]
  nb[,terminal:=node %in% tree$edge[tree$edge[,2]<=ape::Ntip(tree),1]]
  return(nb)
}
