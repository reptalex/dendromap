#' make nodemap for quick referencing of descendants
#' 
#' @export
#' @param tree phylo class object
makeNodeMap <- function(tree){
  D <- length(tree$tip.label)
  nodes <- (D+1):(D+ape::Nnode(tree))
  data.table('node'=nodes,
             'pos'=sapply(nodes,countDescendants,tree),
             'neg'=sapply(nodes,countDescendants,tree,direction='neg'),
             key='node') %>%
    return()
}
