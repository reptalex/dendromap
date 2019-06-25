#' Make basis from phylogenetic tree
#' 
#' @export
#' @param tree \code{phylo} class object
#' @examples
#' library(ape)
#' tree <- rtree(5)
#' treeBasis(tree)
treeBasis <- function(tree){
  getDesc <- function(nd,tr){
    children <- phangorn::Descendants(tr,nd,'children')
    desc <- lapply(children,
                   FUN=function(nd,tr) phangorn::Descendants(tr,nd,'tips')[[1]],
                   tr=tr)
    names(desc) <- paste('node',children,sep='_')
    return(desc)
  }
  basisVec <- function(grp,k){
    v <- rep(0,k)
    r <- length(grp[[1]])
    s <- length(grp[[2]])
    v[grp[[1]]] <- sqrt(s/(r*(r+s)))
    v[grp[[2]]] <- -sqrt(r/(s*(r+s)))
    return(v)
  }
  
  #### find polytomies
  tbl <- table(tree$edge[,1])
  if (tree$Nnode<(length(tree$tip.label)-1)){
    ## polytomies - we omit these nodes as we can't form a phILR basis for them
    polytomies <- as.numeric(names(tbl[tbl>2]))
  } else {
    polytomies <- NULL
  }
  #### find irrelevant nodes
  irrelevant.nodes <- as.numeric(names(tbl[tbl==1]))
  
  #### Obtain set of viable nodes
  k <- ape::Ntip(tree)
  nd <- ape::Nnode(tree)
  nodes <- setdiff(setdiff((k+1):(k+nd),polytomies),irrelevant.nodes)
  if (length(nodes)==0){
    stop('All nodes have either one descendant or more than two descendants')
  }
  
  #### Extract Groups
  Grps <- lapply(nodes,getDesc,tr=tree)
  
  #### make basis
  V <- sapply(Grps,basisVec,k=k)
  colnames(V) <- paste('node',nodes,sep='_')
  rownames(V) <- tree$tip.label
  return(V)
}
