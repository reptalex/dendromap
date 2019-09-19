#' Get index sets of nodes descendant from a focal node, using nodemap
#' @export
#' @param nd node found in \code{nodemap$node}
#' @param nodemap data table produced in \code{\link{makeNodeMap}}
getIndexSets<- function(nd,nodemap){
  x <- vector(mode='list',length=2)
  x[[1]] <- nodemap[node==nd,(nd+1):(nd+pos)*(pos>0)]
  if (all(x[[1]]==0)){
    x[[1]] <- intersect(1,0)
  } 
  
  x[[2]] <- nodemap[node==nd,(nd+pos+1):(nd+pos+neg)*(neg>0)]
  if (all(x[[2]]==0)){
    x[[2]] <- intersect(1,0)
  }
  names(x) <- c('pos','neg')
  return(x)
}