#' Plot treesim object
#' @export
#' @param x treesim class object from \code{\link{treeSim}}
#' @param y optional data matrix. Must have labelled rows containing all tips in \code{x$row.tree} and likewise for \code{x$col.tree}
#' @param color.fcn.clade color function for clade highlights in \code{\link{ggtree}}
#' @param color.fcn.node color function for node circles
#' @param  heatmap.offset offset put into \code{\link{gheatmap}}. Defualt is 0
#' @param col.tr.left Left side boundary of column tree. Default is 0.5
#' @param col.tr.width Width of column tree. Default is 0.45
#' @param col.tr.bottom Bottom of column tree. Deafult is 0.75
#' set.seed(3)
#' m=1e3
#' n=30
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

plot.treesim <- function(x,y=NULL,color.fcn.clade=viridis::viridis,
                         color.fcn.node=viridis::viridis,
                         heatmap.offset=0,
                         col.tr.left=0.5,
                         col.tr.width=0.45,
                         col.tr.bottom=0.75){
  
  vcols <- color.fcn.node(length(unique(S$Paths$col.node)))
  nodecols <- data.table('node'=sort(unique(S$Paths$col.node),decreasing = F),
                         'color'=vcols)
  gtr <- ggtree::ggtree(S$col.tree,branch.length = 'none')+
    ggtree::geom_point2(aes(subset=node %in% nodecols$node),
                        color=nodecols$color,cex=3)+
    ggplot2::coord_flip()+ggplot2::scale_x_reverse()+
    ggplot2::scale_y_reverse()
  
  ###### processing data matrix
  column.order <- gtr$data$label[order(gtr$data$y[1:ape::Ntip(S$col.tree)],decreasing = T)]
  if (is.null(y)){
    X <- S$W %*% S$D %*% t(S$V)
    X <- X[,column.order]
    rownames(X) <- S$row.tree$tip.label
    probs <- binomial(link='logit')$linkinv(X)
    probs[X==0] <- 0.05
    P <- rbinom(nrow(X)*ncol(X),1,probs) %>% 
      matrix(nrow=nrow(X),byrow=F)
    rownames(P) <- rownames(X)
  } else {
    if (!all(S$col.tree$tip.label %in% colnames(y))){
      stop('colnames of y must have all col.tree tip-labels')
    } else if (!all(S$row.tree$tip.label %in% rownames(y))){
      stop('rownames of y must have all row.tree tip-labels')
    } else {
      P <- y[S$row.tree$tip.label,column.order]
    }
  }
  
  gg <- ggtree::ggtree(row.tree,layout='rectangular')
  
  ##flip all the nodes with orientation=-1
  nds <- S$Paths[orientation==-1,row.node]
  for (nd in nds){
    gg <- ggtree::rotate(gg,nd)
  }
  
  ## add clade hilights for basal nodes
  cols <- color.fcn.clade(max(S$Paths$sim))
  ii=0
  start.nodes <- S$Paths[,list(nd=row.node[1]),by=sim]$nd
  for (nd in start.nodes){
    ii=ii+1
    gg <- gg+ggtree::geom_hilight(nd,fill=cols[ii])
  }
  
  Path <- S$Paths
  setkey(Path,row.node)
  gg.cols <- nodecols[match(S$Paths$col.node,node),]$color
  gg <- gg+ggtree::geom_point2(aes(subset=node %in% Path$row.node),
                               color=gg.cols,cex=3)
  
  P <- as.data.frame(P)
  colnames(P) <- colnames(X)
  rownames(P) <- rownames(X)
  gg <- gheatmap(gg,P,color=NA,colnames = FALSE,offset=heatmap.offset)+
    theme(legend.position = 'none')+cowplot::theme_nothing()
  
  output <- cowplot::ggdraw() + 
    cowplot::draw_plot(gtr, x = col.tr.left, y = col.tr.bottom,
                       width = col.tr.width, height = .25)+
    cowplot::draw_plot(gg, x = 0, y = 0.05, width = 1, height = .7)
  
  return(output)
  
}