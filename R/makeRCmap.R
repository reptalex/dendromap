#' make rc-rc mapping for ancestor-descendant relationships & rc tree traversal
#' @export
#' @param rc_table made from \code{\link{makeRCtable}}
#' @param Row_Descendants named \code{getIndexSets} of all row.nodes
#' @param Col_Descendants named \code{getIndexSets} of all col.nodes
makeRCMap <- function(rc_table,Row_Descendants,Col_Descendants){
  row.nodes <- as.numeric(names(Row_Descendants))
  col.nodes <- as.numeric(names(Col_Descendants))
  
  terminalRowNodes <- row.nodes[sapply(Row_Descendants,FUN=function(x) length(unlist(x)))==0]
  terminalColNodes <- row.nodes[sapply(Row_Descendants,FUN=function(x) length(unlist(x)))==0]
  ## we don't have to look for descendants of these row/col nodes
  
  ix <- which((!rc_table$row.node %in% terminalRowNodes) &
                (!rc_table$col.node %in% terminalColNodes))
  maps <- lapply(ix,makeDescendantTable,
                 rc_table,
                 terminalRowNodes,
                 terminalColNodes,
                 Row_Descendants,
                 Col_Descendants)
  RCmap <- rbindlist(maps)
  RCmap[,terminal:=!(descendant %in% ancestor)]
  setkey(RCmap,descendant,ancestor)
  return(RCmap)
}
