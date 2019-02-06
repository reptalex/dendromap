#' Make basis from phylogenetic tree
#' 
#' @export
#' @param tree \code{phylo} class object
#' @examples
#' library(ape)
#' tree <- rtree(5)
#' treeBasis(tree)
treeBasis <- function(tree){
  k <- ape::Ntip(tree)
  nd <- ape::Nnode(tree)
  getDesc <- function(nd,tr){
    children <- phangorn::Descendants(tr,nd,'children')
    desc <- lapply(children,
                   FUN=function(nd,tr) phangorn::Descendants(tr,nd,'tips')[[1]],
                   tr=tr)
    names(desc) <- paste('node',children,sep='_')
    return(desc)
  }
  Grps <- lapply((k+1):(k+nd),getDesc,tr=tree)
  
  basisVec <- function(grp,k){
    v <- rep(0,k)
    r <- length(grp[[1]])
    s <- length(grp[[2]])
    v[grp[[1]]] <- sqrt(s/(r*(r+s)))
    v[grp[[2]]] <- -sqrt(r/(s*(r+s)))
    return(v)
  }
  
  V <- sapply(Grps,basisVec,k=k)
  colnames(V) <- paste('node',(k+1):(k+nd),sep='_')
  rownames(V) <- tree$tip.label
  return(V)
}
