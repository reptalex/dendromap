#' Internal functin helping makeRCmap find rc-rc ancestor-descendant relations
#' @export
#' @param i index of rc_table whose descendants are found
#' @param terminalRowNodes row nodes in rc_table with no descendants
#' @param terminalColNodes same for col nodes in col tree
#' @param Row_Descendants list of descendants for each row node
#' @param Col_Descendants list of descendants for each col node
makeDescendantTable <- function(i=1,rc_table,
                                terminalRowNodes,
                                terminalColNodes,
                                Row_Descendants,
                                Col_Descendants){
  row.node <- rc_table$row.node[i]
  col.node <- rc_table$col.node[i]
  if (row.node %in% terminalRowNodes | 
      col.node %in% terminalColNodes){
    return(NULL)
  } else {
    sgn <- sign(rc_table$stat[i])
    if (sgn==1){ ## pos-pos and neg-neg
      rowdescPos <- Row_Descendants[[toString(row.node)]][['pos']]
      coldescPos <- Col_Descendants[[toString(col.node)]][['pos']]
      rowdescNeg <- Row_Descendants[[toString(row.node)]][['neg']]
      coldescNeg <- Col_Descendants[[toString(col.node)]][['neg']]
    } else {
      rowdescPos <- Row_Descendants[[toString(row.node)]][['pos']]
      coldescPos <- Col_Descendants[[toString(col.node)]][['neg']]
      rowdescNeg <- Row_Descendants[[toString(row.node)]][['neg']]
      coldescNeg <- Col_Descendants[[toString(col.node)]][['pos']]
    }
    ix <- rc_table[(row.node %in% rowdescPos & col.node %in% coldescPos) |
                     (row.node %in% rowdescNeg & col.node %in% coldescNeg)]$rc_index
    if (length(ix)>0){
      return(data.table('ancestor'=rc_table$rc_index[i],
                        'descendant'=ix))
    } else {
      return(NULL)
    }
  }
}