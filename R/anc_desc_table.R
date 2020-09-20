#' Get ancestor-descendant relations a given rc pair
#' @export
#' @param i integer: row of \code{rc_table} for extraction of ancestor and descendant \code{rc_index} relations
#' @param rc_table see \code{\link{make_rc_table}}
#' @param rowTerminals Indexes of terminal row tree edges
#' @param colTerminals Indexes of terminal column tree edges
#' @param rowDescendants List of descendants made from \code{edge_registry}. Element \code{i} contains all edges descendant from edge \code{i} in the row tree
#' @param colDescendants List of descendants made from \code{edge_registry}. Element \code{j} contains all edges descendant from edge \code{j} in the column tree
#' @examples
#' 
#' # none yet - for internal dendromap use
anc_desc_table <- function(i=1,rc_table,rowTerminals,colTerminals,
                         rowDescendants,colDescendants){
    row.edge <- rc_table$row.edge[i]
    col.edge <- rc_table$col.edge[i]
    if (row.edge %in% rowTerminals|
        col.edge %in% colTerminals){
      return(NULL)
    } else {
      
      rowDesc <- rowDescendants[[row.edge]]
      colDesc <- colDescendants[[col.edge]]
      
      ixs <- rc_table[(row.edge %in% rowDesc) & (col.edge %in% colDesc)]$rc_index
      if (length(ixs)>0){
        return(data.table('ancestor'=rc_table[i]$rc_index,
                          'descendant'=ixs))
      } else {
        return(NULL)
      }
    }
}