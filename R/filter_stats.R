#' Plot dendromap object
#' @export
#' @param x output from \code{\link{lineage_stats}}
#' @param Lineages see \code{\link{get_lineages_and_stats}}
#' @param rc_table see \code{\link{make_rc_table}}
#' @param rc_relations see \code{\link{get_rc_relations}}
#' @param rowDescendants List of descendants made from \code{edge_registry}. Element \code{i} contains all edges descendant from edge \code{i} in the row tree
#' @param row.tree phylo class object
#' @examples
filter_stats = function(x,Lineages,rc_table,rc_relations,rowDescendants,row.tree){
  x <- x[order(F_stat,decreasing = T)]
  
  i=1
  while(i<=nrow(x)){
    id <- x$lineage_id[i]
    row.edg <- rc_table[rc_index==id,row.edge]
    descendants <- rowDescendants[[row.edg]]
    ancestors <- edge_ancestors(row.edg,row.tree)
    incompatible_lineages <- setdiff(Lineages[row.edge %in% c(descendants,ancestors),lineage_id],id)
    x <- x[!lineage_id %in% incompatible_lineages]
    i=i+1
  }
  return(x)
}
