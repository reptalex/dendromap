#' label children of a node to keep track of tree orientation
#' @export
#' @param node internal node of \code{tree}
#' @param tree \code{phylo} class object
#' @examples
#' set.seed(1)
#' library(ape)
#' tree <- rtree(10)
#' labelledChildren(11,tree)
labelledChildren <- function(node,tree){
  ch <- phangorn::Descendants(tree,node,'children')
  names(ch) <- c('a','b')
  return(ch)
}
