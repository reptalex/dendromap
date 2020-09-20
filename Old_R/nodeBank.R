#' Make a dataset of nodes in a tree and whether/not they are manifested
#' 
#' @export
#' @param node node in tree
#' @param tree \code{phylo} class object
#' @param prob probability between 0 and 1
#' @examples
#' set.seed(1)
#' library(ape)
#' tree <- rtree(10)
#' nb <- nodeBank(11,tree,0.5)
nodeBank <- function(node,tree,prob=1,propensity=NULL,include.node=FALSE){
  nb <- data.table('node'=phangorn::Descendants(tree,node,'all'))
  if (include.node){
    nb <- rbind(data.table('node'=node),nb)
  }
  nb <- nb[node %in% (ape::Ntip(tree)+1:tree$Nnode)]
  if (is.null(propensity)){
    nb[,propensity:=-log(1-prob)]
  } else {
    nb[,propensity:=propensity]
    prob <- 1-exp(-propensity)
  }
  nb[,manifest:=sample(c(T,F),size=.N,replace = T,prob=c(prob,1-prob))]
  # tree$edge[tree$edge[,2]<=ape::Ntip(tree),1]
  terminal.nodes <- table(tree$edge[tree$edge[,2]<=ape::Ntip(tree),1])
  terminal.nodes <- as.numeric(names(which(terminal.nodes==2)))
  nds <- data.table('node'=1:(ape::Ntip(tree)+ape::Nnode(tree)),
                    'depth'=ape::node.depth.edgelength(tree))
  nds <- nds[node %in% nb$node]
  setkey(nds,node)
  setkey(nb,node)
  nb <- nb[nds]
  nb[,terminal:=node %in% terminal.nodes]
  
  return(nb)
}
