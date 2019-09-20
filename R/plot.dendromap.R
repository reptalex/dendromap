#' Plot dendromap object
#' @export
#' @param x dendromap class object from \code{\link{treeSim}} or \code{\link{dendromap}}
#' @param y optional data matrix. Must have labelled rows containing all tips in \code{x$row.tree} and likewise for \code{x$col.tree}
#' @param color.fcn.clade color function for clade highlights in \code{\link{ggtree}}
#' @param color.fcn.node color function for node circles
#' @param  heatmap.offset offset put into \code{\link{gheatmap}}. Defualt is 0
#' @param col.tr.left Left side boundary of column tree. Default is 0.5
#' @param col.tr.width Width of column tree. Default is 0.45
#' @param col.tr.bottom Bottom of column tree. Deafult is 0.75
#' @examples
#' set.seed(3)
#' m=1e3
#' n=7
#' row.tree <- rtree(m) %>% phytools::force.ultrametric()
#' col.tree <- rtree(n)
#' 
#' S <- treeSim(10,row.tree,col.tree,prob.row=0.7,prob.col=0.8,
#'              col.node = n+1,fix.col.node = T,sd = 1e3,
#'              row.depth.min=2,row.depth.max=3)
#' 
#' plot(S,col.tr.left = 0.47,
#'      col.tr.width = 0.505,
#'      col.tr.bottom = 0.74)
plot.dendromap <- function(x,y=NULL,color.fcn.clade=viridis::viridis,
                         color.fcn.node=viridis::viridis,
                         heatmap.offset=0,
                         col.tr.left=0.5,
                         col.tr.width=0.45,
                         col.tr.bottom=0.75){
  
  if (is.null(x$Lineages$orientation)){
    x$Lineages[,orientation:=sign(stat)]
  }
  vcols <- color.fcn.node(length(unique(x$Lineages$col.node)))
  nodecols <- data.table('node'=sort(unique(x$Lineages$col.node),decreasing = F),
                         'color'=vcols)
  gtr <- ggtree::ggtree(x$col.tree,branch.length = 'none')+
    ggtree::geom_point2(ggplot2::aes(subset=node %in% nodecols$node),
                        color=nodecols$color,cex=3)+
    ggplot2::coord_flip()+ggplot2::scale_x_reverse()+
    ggplot2::scale_y_reverse()
  
  ###### processing data matrix
  column.order <- gtr$data$label[order(gtr$data$y[1:ape::Ntip(x$col.tree)],decreasing = T)]
  if (is.null(y)){
    if (is.null(x$Data)){
      if (is.null(x$W) | is.null(x$V) | is.null(x$D)){
        y <- base::matrix(0,nrow=ape::ntip(x$row.tree),ncol=ape::ntip(col.tree))
      } else {
        y <- x$W %*% x$D %*% t(x$V)
      }
    } else {
      y <- x$Data
    }
  }
  
  gg <- ggtree::ggtree(x$row.tree,layout='rectangular',branch.length = 'none')
  
  ##flip all the nodes with orientation=-1
  nds <- x$Lineages[orientation==-1,row.node]
  for (nd in nds){
    gg <- ggtree::rotate(gg,nd)
  }
  
  ## add clade hilights for basal nodes
  cols <- color.fcn.clade(max(x$Lineages$Lineage))
  ii=0
  start.nodes <- x$Lineages[,list(nd=min(row.node)),by=Lineage]$nd
  for (nd in start.nodes){
    ii=ii+1
    gg <- gg+ggtree::geom_hilight(nd,fill=cols[ii])
  }
  
  Path <- x$Lineages
  setkey(Path,row.node)
  gg.cols <- nodecols[match(x$Lineages$col.node,node),]$color
  gg <- gg+ggtree::geom_point2(ggplot2::aes(subset=node %in% Path$row.node),
                               color=gg.cols,cex=3)
  
  gg <- ggtree::gheatmap(gg,as.matrix(y),color=NA,colnames = FALSE,offset=heatmap.offset)+
    ggtree::theme(legend.position = 'none')+cowplot::theme_nothing()
  
  output <- cowplot::ggdraw() + 
    cowplot::draw_plot(gtr, x = col.tr.left, y = col.tr.bottom,
                       width = col.tr.width, height = .25)+
    cowplot::draw_plot(gg, x = 0, y = 0.05, width = 1, height = .7)
  
  return(output)
}
