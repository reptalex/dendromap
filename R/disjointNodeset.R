#' Find nodes not on rootpaths between a set of nodes and their descendants
#' @export
#' @param nodes nodes of input \code{tree}
#' @param tree \code{phylo} class object
#' @examples
#' set.seed(1)
#' library(ape)
#' tree <- rtree(10)
#' nodes <- c(15,19)
#' other.nodes <- disjointNodeset(nodes,tree)
#' plot(tree)
#' nodes <- setdiff(Ntip(tree)+1:tree$Nnode,other.nodes)
#' nodelabels(nodes,nodes,bg=rgb(0,.5,.5,1))
#' nodelabels(c(15,19),c(15,19),col='red')
#' nodelabels(other.nodes,other.nodes,bg=rgb(.9,.9,0,1))

disjointNodeset <- function(nodes,tree){
  nds <- phangorn::Descendants(tree,nodes,'all') %>%
    unlist %>% unique
  rootpaths <- sapply(nds,FUN=function(nd,tr) nodepath(tr,nd,ape::Ntip(tr)+1),tree) %>%
    unlist %>% unique
  nds <- unique(c(nds,rootpaths))
  setdiff(1:ape::Nnode(tree)+ape::Ntip(tree),nds) %>%
    return()
}