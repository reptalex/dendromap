#' Build lineage from a given rc_index, using only P <= p_thresh
#' @export
#' @param basal_ix index found in \code{rc_table$rc_index}
#' @param p_thresh Positive numeric, threshold P-value for inclusion of rc-pairs
#' @param rc_relations 
#' @param rc_table see \code{\link{make_rc_table}}
#' @param row.tree \code{phylo} class object
#' @param rowDescendants List of descendants made from \code{edge_registry}. Element \code{i} contains all edges descendant from edge \code{i} in the row tree
#' @param colDescendants List of descendants made from \code{edge_registry}. Element \code{j} contains all edges descendant from edge \code{j} in the column tree
#' @examples
get_lineage <- function(basal_ix,p_thresh=1,rc_relations.=rc_relations,rc_table.=rc_table,row.tree.=row.tree,
                        rowDescendants.=rowDescendants,colDescendants.=colDescendants,...){
  rct <- rc_table[P<=p_thresh]
  rcm <- rc_relations[max_P<=p_thresh]
  basal_row_edge <- rct[rc_index==basal_ix,row.edge]
  basal_col_edge <- rct[rc_index==basal_ix,col.edge]
  row_descs <- rowDescendants[[basal_row_edge]]
  col_descs <- colDescendants[[basal_col_edge]]
  rc_tbl <- rct[row.edge %in% c(basal_row_edge,row_descs) & 
                  col.edge %in% c(basal_col_edge,col_descs)]
  lineage <- assemble_lineage(rc_tbl,rcm,rt,row.tree)
  lineage[,lineage_id:=basal_ix]
  setkey(lineage,col.edge,P)
  return(lineage)
}
