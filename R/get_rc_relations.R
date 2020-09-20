#' Get all ancester-descendant relations of rc pairs
#' @export
#' @param rc_table see \code{\link{make_rc_table}}
#' @param rowDescendants List of descendants made from \code{edge_registry}. Element \code{i} contains all edges descendant from edge \code{i} in the row tree
#' @param colDescendants List of descendants made from \code{edge_registry}. Element \code{j} contains all edges descendant from edge \code{j} in the column tree
#' @examples
get_rc_relations <- function(rc_table,rowDescendants,colDescendants){
  
  rowTerminals <- which(sapply(rowDescendants,FUN=function(x) length(x)==0))
  colTerminals <- which(sapply(colDescendants,FUN=function(x) length(x)==0))
  
  ix <- which((!rc_table$row.edge %in% rowTerminals) &
                (!rc_table$col.edge %in% colTerminals))
  maps <- lapply(ix,anc_desc_table,
                 rc_table,
                 rowTerminals,
                 colTerminals,
                 rowDescendants,
                 colDescendants)
  rc_ancestry_relations <- rbindlist(maps)
  rc_ancestry_relations[,terminal:=!(descendant %in% ancestor)]
  setkey(rc_ancestry_relations,descendant,ancestor)
  return(rc_ancestry_relations)
}