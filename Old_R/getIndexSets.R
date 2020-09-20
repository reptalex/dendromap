#' Get index sets of nodes descendant from a focal node, using map
#' @export
#' @param x node found in \code{map$node}
#' @param map data table produced in \code{\link{makemap}}
getIndexSets<- function(x,map,method='node'){
  if (method=='node'){
    y <- vector(mode='list',length=2)
    y[[1]] <- map[node==x,(x+1):(x+pos)*(pos>0)]
    if (all(y[[1]]==0)){
      y[[1]] <- intersect(1,0)
    } 
    
    y[[2]] <- map[node==x,(x+pos+1):(x+pos+neg)*(neg>0)]
    if (all(y[[2]]==0)){
      y[[2]] <- intersect(1,0)
    }
    names(y) <- c('pos','neg')
  } else {
    y <- setdiff(map[edge==x,seq(x,x+n)],x)
  }
  return(y)
}
