#' Visualize a single, manifested path on two trees
#' @export
#' @param Path output from \code{\link{treeMap}}
#' @param col.tree \code{phylo} class object from \code{treeMap}
#' @param row.tree \code{phylo} class object from \code{treeMap}
#' @examples
#' set.seed(1)
#' library(ape)
#' row.tree <- rtree(100)
#' col.tree <- rtree(10)
#' S <- treeMap(row.tree,col.tree,101,11,0.2,1)
#' pathViz(S,col.tree,row.tree)
pathViz <- function(Path,col.tree,row.tree){
  colmap <- viridis::viridis(nrow(Path))
  par(mfrow=c(1,2))
  plot(row.tree,main='Row tree')
  nodelabels(Path$row.node,Path$row.node,bg = colmap,cex=2,frame = 'circle')
  plot(col.tree,main='Column tree')
  nodelabels(Path$col.node,Path$col.node,bg = colmap,cex=2,frame='circle')
}