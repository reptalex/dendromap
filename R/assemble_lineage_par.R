#' Plot dendromap object
#' @export
#' @param rc_tbl see \code{\link{make_rc_table}}. Subset of \code{rc_table}
#' @param rc_relations see \code{\link{get_rc_relations}}
#' @param rc_table see \code{\link{make_rc_table}}
#' @param row.tree \code{phylo} class object
#' @param rowDescendants List of descendants made from \code{edge_registry}. Element \code{i} contains all edges descendant from edge \code{i} in the row tree
#' @param colDescendants List of descendants made from \code{edge_registry}. Element \code{j} contains all edges descendant from edge \code{j} in the column tree
#' @examples
assemble_lineage_par <- function(rc_tbl,rc_relations.=rc_relations,rc_table.=rc_table,row.tree.=row.tree,
                             rowDescendants.=rowDescendants,colDescendants.=colDescendants){
  ### RULES:
  ## 1) remove col.edges with more significant values among descendants (implies signal lost with inclusion of sister taxa)
  ##     see function basal_dominance_filter - this is done externally
  ## 2) pick most basal row.edge, choosing col.edge with lowest P-value  -- from filtering, this is guaranteed to be greater than that of any descendant for same col.edge
  ## 3) remove row.edge and all descendant rc_indexes with same col.edge from rc_tbl
  ## 4) repeat 2-3 until rc_tbl is empty
  
  lineage=NULL
  done=F
  n=0
  setkey(rc_tbl,row.edge,P)
  while(!done){
    n=n+1
    lineage <- rbind(lineage,rc_tbl[1,]) ## next most basal row edge, and its lowest P-value col.edge
    ix <- rc_tbl$rc_index[1]
    incompatibles <- c(ix,incompatible_descendants(ix,rc_tbl,rowDescendants,colDescendants))
    rc_tbl <- rc_tbl[!rc_index %in% incompatibles]
    if (nrow(rc_tbl)==0){
      done=TRUE
    }
  }
  lineage <- merge_twin_sisters(lineage,row.tree,rc_table)
  return(lineage)
}