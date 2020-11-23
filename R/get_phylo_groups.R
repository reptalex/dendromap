#' Get list of groups defined by phylogeny
#' @export
#' @param tree Phylogeny
#' @return list of groups and their complements, i.e. lists of groups split by each edge of the unrooted phylogeny
#' @useDynLib dendromap
#' @useDynLib dendromap bifurcations
#' @examples
#' tr <- rtree(8)
#' tr$tip.label <- 1:8
#' par(mfrow=c(1,1))
#' plot.phylo(tr,use.edge.length=FALSE);edgelabels()
#' Grps <- get_phylo_groups(tr)

######################### getPhyloGroups ######################################

get_phylo_groups <- function(tree) {
  .fillInRes <- function(i) {
    res <- vector('list', length=2)
    res[[1]] <- which(edge.tips[i, ] == 1)
    res[[2]] <- which(edge.tips[i, ] == 0)
    res
  }
  edge.tips <- .Call("bifurcations", PACKAGE="dendromap",
                     as.integer(tree$edge[, 1]),
                     as.integer(tree$edge[, 2]),
                     as.integer(length(tree$tip.label)))
  lapply(1:nrow(edge.tips), .fillInRes)
}
