#' Visualize a single, manifested path on two trees
#' @export
#' @param Path output from \code{\link{treeMap}}
#' @param row.tree \code{phylo} class object from \code{treeMap}
#' @param col.tree \code{phylo} class object from \code{treeMap}
#' @param ... additional input args to plot and nodelabels
#' @examples
#' set.seed(1)
#' library(ape)
#' row.tree <- rtree(100)
#' col.tree <- rtree(10)
#' S <- treeMap(row.tree,col.tree,101,11,0.2,1)
#' pathViz(S,row.tree,col.tree)
pathViz <- function(Path,row.tree,col.tree,...){
  
  colmap <- data.table('col.node'=sort(unique(Path$col.node),decreasing = F),
                       'col'=viridis::viridis(length(unique(Path$col.node))))
  
  par(mfrow=c(1,2))
  plot(row.tree,main='Row tree',...)
  nodelabels(Path$row.node,Path$row.node,
             bg = colmap[match(Path$col.node,colmap$col.node),col],
             frame = 'circle',...)
  
  plot(col.tree,main='Column tree',...)
  nodelabels(colmap$col.node,colmap$col.node,bg = colmap$col,frame='circle',...)
}
