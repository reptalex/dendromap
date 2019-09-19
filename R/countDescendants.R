#' Count descendants of a node
#' 
#' @export
#' @param node numeric node of \code{tree} which has only two daughter clades
#' @param tree phylo class object
#' @param direction characer, either "pos" or "neg", indicating the direction of descendants to count (e.g. which daughter clade to count)/
countDescendants <- function(node,tree,direction='pos'){
  children <- phangorn::Descendants(tree,node,type='children')
  if (length(children)<2){
    n=0
  } else {
    D <- length(tree$tip.label)
    if (direction=='pos'){
      if (children[1]<=D){
        n <- 0
      } else {
        n <- length(setdiff(phangorn::Descendants(tree,children[1],type = 'all'),1:D))+1
      }
    } else {
      if (children[2] <=D){
        n <- 0
      } else {
        n <- length(setdiff(phangorn::Descendants(tree,children[2],type = 'all'),1:D))+1
      }
    }
  }
  return(n)
}
