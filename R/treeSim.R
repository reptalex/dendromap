#' Simulate multiple treeMaps on row.tree
#' @export
#' @param n max number of independent lineages on \code{row.tree} to map to column tree
#' @param row.tree \code{phylo} class object
#' @param col.tree \code{phylo} class object
#' @param prob.row probability of a node in the row.tree manifesting a map to col.tree. Default is 1
#' @param prob.col probability of a node in the col.tree manifesting a map to row.tree. Default is 1
#' @param sd optional positive numeric, standard deviation input to \code{\link{rnorm}} to create output diagonal matrix
#' @param W optional matrix, \code{\link{treeBasis}} of \code{row.tree}
#' @param V optional matrix, \code{\link{treeBasis}} of \code{col.tree}
#' @param row.node optional start node for row.tree
#' @param col.node optional start node for col.tree
#' @param row.nodeset optional set of nodes to consider for random start of \code{treeMap}
#' @param col.nodeset optional set of nodes to consider for random start of \code{treeMap}
#' @param row.depth.min optional numeric minimum depth for nodes for random start of \code{treeMap}
#' @param row.depth.max optional numeric maximum depth for nodes for random start of \code{treeMap}
#' @param col.depth.min optional numeric minimum depth for nodes for random start of \code{treeMap}
#' @param col.depth.max optional numeric maximum depth for nodes for random start of \code{treeMap}
#' @param fix.col.node logical: whether or not to fix the column node for the start of every \code{treeMap}
#' @examples 
#' set.seed(3)
#' m=1e3
#' n=30
#' row.tree <- rtree(m) %>% phytools::force.ultrametric()
#' col.tree <- rtree(n)
#' 
#' S <- treeSim(10,row.tree,col.tree,prob.row=0.7,prob.col=0.8,
#'              col.node = n+1,fix.col.node = T,sd = 1e3,
#'              row.depth.min=2,row.depth.max=3)
#' S <- treeSim(10,row.tree,col.tree,prob.row=0.7,prob.col=0.8,sd = 1e3,
#'              row.depth.min=2,row.depth.max=3)              
#' 
#' plot(S,col.tr.left = 0.47,
#'      col.tr.width = 0.505,
#'      col.tr.bottom = 0.74)

treeSim <- function(n,row.tree,col.tree,
                    prob.row=1,prob.col=1,sd=1,W=NULL,V=NULL,
                    row.node=NULL,col.node=NULL,
                    row.nodeset=NULL,col.nodeset=NULL,
                    row.depth.min=NULL,row.depth.max=NULL,
                    col.depth.min=NULL,col.depth.max=NULL,
                    fix.col.node=FALSE){
  
  if (fix.col.node==TRUE & is.null(col.node)){
    col.node <- 1+ape::Ntip(col.tree)
  }
  
  if (is.null(row.depth.min)){
    row.depth.min <- 0
  }
  if (is.null(row.depth.max)){
    row.depth.max <- Inf
  }
  if (is.null(col.depth.min)){
    col.depth.min <- 0
  }
  if (is.null(col.depth.max)){
    col.depth.max <- Inf
  }
  
  depths <- ape::node.depth.edgelength(row.tree)
  if (is.null(row.nodeset)){
    row.nodeset <- intersect(1:ape::Nnode(row.tree)+length(row.tree$tip.label),
                             which(depths<=row.depth.max & depths>=row.depth.min))
  } else {
    row.nodeset <- intersect(row.nodeset,
                             which(depths<=row.depth.max & depths>=row.depth.min))
  }
  
  if (is.null(col.nodeset)){
    col.nodeset <- intersect(1:ape::Nnode(col.tree)+length(col.tree$tip.label),
                             which(depths<=col.depth.max & depths>=col.depth.min))
  } else {
    col.nodeset <- intersect(col.nodeset,
                             which(depths<=col.depth.max & depths>=col.depth.min))
  }
  
  
  Path <- vector(mode='list',length=n)
  for (ii in 1:n){
    if (ii==1){
      if (is.null(row.node)){
        row.node <- sample(row.nodeset,1)
      }
      if (is.null(col.node)){
        col.node <- sample(col.nodeset,1)
      }
      
      Path[[ii]] <- treeMap(row.tree,col.tree,row.node,col.node,prob.row,prob.col)
    } else {
      row.start <- Path[[ii-1]]$row.node[1]
      row.nodeset <- intersect(row.nodeset,disjointNodeset(row.start,row.tree))
      
      if (length(row.nodeset)>0){
        if (!fix.col.node){
          col.node <- sample(col.nodeset,1)
        }
        row.node <- sample(row.nodeset,1)
        Path[[ii]] <-  treeMap(row.tree,col.tree,row.node,col.node,prob.row,prob.col)
      } else {
        Path <- Path[1:(ii-1)]
        break
      }
    }
  }
  for (i in 1:length(Path)){
    Path[[i]][,sim:=i]
  }
  output <- NULL
  output$Paths <- data.table::rbindlist(Path)
  if (is.null(W)){
    W <- treeBasis(row.tree)
  }
  if (is.null(V)){
    V <- treeBasis(col.tree)
  }
  
  output$W <- W[,output$Paths$row.node-ape::Ntip(row.tree)]
  dd <- rnorm(nrow(output$Paths),sd=sd)
  dd <- abs(dd)*output$Paths$orientation
  output$D <- diag(dd)
  output$V <- V[,output$Paths$col.node-ape::Ntip(col.tree)]
  output$row.tree <- row.tree
  output$col.tree <- col.tree
  class(output) <- 'treesim'
  return(output)
}
