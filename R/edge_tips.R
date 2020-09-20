#' Table of edges and max/min values of descendant tips
#' @export
#' @param tree \code{phylo} class object
#' @examples
edge_tips <- function(tree){
  grps <- phylofactor::getPhyloGroups(tree)
  mat <- lapply(grps,FUN=function(gg)
    data.table('min'=min(gg[[1]]),'max'=max(gg[[1]])))
  mat <- rbindlist(mat)
  mat[,edge:=1:.N]
  setkey(mat,edge)
  mat[,tip:=(min==max)]
  return(mat)
}