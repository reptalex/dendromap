#' make rc-rc mapping for ancestor-descendant relationships & rc tree traversal
#' @export
#' @param rc_table made from \code{\link{makeRCtable}}
#' @param row.nodemap made from \code{\link{makeNodeMap}} of row.tree
#' @param col.nodemap made from \code{\link{makeNodeMap}} of col.tree
makeRCMap <- function(rc_table,row.nodemap,col.nodemap){
  row.nodes <- unique(rc_table$row.node)
  col.nodes <- unique(rc_table$col.node)
  Row_Descendants <- lapply(row.nodes,getIndexSets,row.nodemap) %>%
    lapply(FUN=function(x,a) lapply(x,intersect,a),a=row.nodes)
  Col_Descendants <- lapply(col.nodes,getIndexSets,col.nodemap) %>%
    lapply(FUN=function(x,a) lapply(x,intersect,a),a=col.nodes)
  
  names(Row_Descendants) <- row.nodes
  names(Col_Descendants) <- col.nodes
  
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