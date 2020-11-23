#' Create edge contrast vectors for tree
#' @export
#' @param tree phylo class object
#' @example
#' tr <- rtree(5)
#' edge_contrasts(tr)
edge_contrasts <- function(tree){
  get_phylo_groups(tree) %>% sapply(contrast_vector,n=length(tree$tip.label))
}