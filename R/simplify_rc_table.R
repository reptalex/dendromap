#' Remove (r,c) pairs without ancestors/descendants, and impose basal dominance
#' @export
#' @param rc_table see \code{\link{make_rc_table}}
#' @param rc_relations see \code{\link{get_rc_relations}}
#' @param rowDescendants List of descendants made from \code{edge_registry}. Element \code{i} contains all edges descendant from edge \code{i} in the row tree
#' @param colDescendants List of descendants made from \code{edge_registry}. Element \code{j} contains all edges descendant from edge \code{j} in the column tree
#' @examples
simplify_rc_table <- function(rc_table,rc_relations,rowDescendants,colDescendants){
  ancs <- unique(rc_relations$ancestor)
  descs <- unique(rc_relations$descendant)
  basal_indexes <- setdiff(ancs,descs) ## ancestors w/o descendant
  
  ### We'll filter rc_table to: (1) only edges with ancestor/descendant relation
  ### i.e. must be in rc_relations
  rc_table <- rc_table[rc_index %in% c(ancs,descs)]
  setkey(rc_table,row.edge,P)
  rc_table <- basal_dominance_filter(rc_tbl=rc_table,rowDescendants)
  return(rc_table)
}