#' Warps tree
#' @export
#' @param tree phylo class object
#' @param trans character, either 'log' or 'power', indicating whether to use \code{a*log(x/a+1)} or \code{x^a/a} warping
#' @param alpha see functional forms for \code{trans}. Default for 'log' is 1, default for 'power' is 1/2
#' @examples 
#' set.seed(1)
#' tr <- rtree(100)
#' tr <- phytools::force.ultrametric(tr)
#' par(mfrow=c(2,1))
#' plot(tr,show.tip.label = F)
#' plot(warp.phylo(tr),show.tip.label = F)
warpTree = function(tree,trans='log',alpha=NULL){
  edges <- tree$edge
  nodes_depth <- node.depth.edgelength(tree)
  nodes_depth <- max(nodes_depth)-nodes_depth  ##invert height to depth so root is deepest
  
  if (trans=='log'){
    if (is.null(alpha)){
      alpha=1
    }
    y <- function(x,alpha.=alpha) alpha*log(x/alpha+1)
  } else {
    if (is.null(alpha)){
      alpha=1/2
    }
    y <- function(x,alpha.=alpha) x^alpha/alpha
  }
  
  for (i in 1:nrow(edges)){
    beg <- nodes_depth[edges[i, 1]]
    end <- nodes_depth[edges[i, 2]]
    tree$edge.length[i] <- y(beg)-y(end)
  }
  return(tree)
}
