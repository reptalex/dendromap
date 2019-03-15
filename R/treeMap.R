#' Map events on column tree to nodes on row tree
#' @export
#' @param row.tree \code{phylo} class object
#' @param col.tree \code{phylo} class object
#' @param row.node integer. optional node for starting treeMap. Must be a node in \code{row.tree}
#' @param col.node integer. optional node for starting treeMap. Must be a node in \code{col.tree}
#' @param prob.row probability, between 0 and 1, that a descendant node is manifested in the tree map.
#' @param prob.col same as \code{prob.row} but for nodes in column tree.
#' @param col.nb optional \code{\link{nodeBank}} with defined propensities for column node manifestation
#' @examples
#' set.seed(1)
#' library(ape)
#' library(twotree)
#' library(ggplot2)
#' row.tree <- rtree(100)
#' col.tree <- rtree(10)
#' S <- treeMap(row.tree,col.tree,101,11,0.2,1)
#' pathViz(S,row.tree,col.tree)
#' 
#' S2 <- treeMap(row.tree,col.tree,101,11,0.2,1,use.depths=T)
#' pathViz(S2,row.tree,col.tree)
#' 
#' ### We can also simulate a particular column node having a
#' ### higher propensity:
#' col.nb <- nodeBank(11,col.tree,propensity=2)
#' col.nb[node %in% c(13,17,18),propensity:=Inf]  ## high chance of these nodes
#' 
#' S3 <- treeMap(row.tree,col.tree,row.node=101,col.node=11,prob.row=0.1,
#'               col.nb=col.nb,use.depths=T)
#' pathViz(S3,row.tree,col.tree)
#' 
#' set.seed(1)
#' library(ape)
#' row.tree <- rtree(100)
#' col.tree <- rtree(10)
#' 
#' set.seed(1)
#' library(ape)
#' row.tree <- rtree(1000)
#' col.tree <- rtree(300)
#' 
#' S <- treeMap(row.tree,col.tree,row.node=1001,col.node=301,prob.row=0.1,
#'               col.nb=nodeBank(301,col.tree,propensity=1),use.depths=T)
#' row.depths <- data.table('row.depth'=node.depth.edgelength(row.tree))
#' row.depths[,row.node:=1:.N]
#' setkey(row.depths,row.node)
#' col.depths <- data.table('col.depth'=node.depth.edgelength(col.tree))
#' col.depths[,col.node:=1:.N]
#' setkey(col.depths,col.node)
#' setkey(S,row.node)
#' S <- row.depths[S]
#' setkey(S,col.node)
#' S <- col.depths[S]
#' 
#' ggplot(S,aes(row.depth,col.depth))+geom_point()+geom_abline(intercept=0,slope=1)
treeMap <- function(row.tree,col.tree,row.node=NULL,col.node=NULL,
                    prob.row=1,prob.col=1,col.nb=NULL,use.depths=F){
  if (is.null(row.node)){
    row.node <- ape::Ntip(row.tree)+1
  }
  if (is.null(col.node)){
    col.node <- ape::Ntip(col.tree)+1
  }
  
  
  row.nb <- nodeBank(ape::Ntip(row.tree)+1,row.tree,prob.row)
  if (is.null(col.nb)){
    col.nb <- nodeBank(ape::Ntip(col.tree)+1,col.tree,prob.col)
    col.nb[,propensity:=1]
  } else {
    if (is.null(col.nb$propensity)){
      col.nb[,propensity:=1]
    }
  }
  Path <- data.table('row.node'=row.node,
                     'col.node'=col.node,
                     'orientation'=sign(rnorm(1)),
                     'terminated'=nrow(row.nb)==0 | nrow(col.nb)==0)

  while (any(!Path$terminated)){
    ix <- min(which(!Path$terminated))
    ch.row <- labelledChildren(Path$row.node[ix],row.tree)
    ch.col <- labelledChildren(Path$col.node[ix],col.tree)
    orientation <- Path$orientation[ix]
    if (length(ch.col)>0 & length(ch.row)>0){
      
      colA <- manifestChildren(ch.col['a'],col.tree,col.nb,use.depths)
      colB <- manifestChildren(ch.col['b'],col.tree,col.nb,use.depths)
      
      if (length(colA)>0){
        if (orientation>0){
          Path <- rbind(Path,manifestRowChildren(ch.row['a'],row.tree,row.nb,col.nb,colA,use.depths))
        } else {
          Path <- rbind(Path,manifestRowChildren(ch.row['b'],row.tree,row.nb,col.nb,colA,use.depths))
        }
      }
      if (length(colB)>0){
        if (orientation>0){
          Path <- rbind(Path,manifestRowChildren(ch.row['b'],row.tree,row.nb,col.nb,colB,use.depths))
        } else {
          Path <- rbind(Path,manifestRowChildren(ch.row['a'],row.tree,row.nb,col.nb,colB,use.depths))
        }
      }
      Path$terminated[ix] <- TRUE
    } else {
      Path$terminated[ix] <- TRUE
    }
  }
  return(Path)
}
