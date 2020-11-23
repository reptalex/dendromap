#' Find descendant rc_indexes incompatible with input rc_index
#' @export
#' @param rc_ix index found in \code{rc_tbl$rc_index}
#' @param rc_tbl see \code{\link{make_rc_table}}. Internally in \code{\link{dendromap}}, this is a trimmed version of the global \code{rc_table}
#' @param rowDescendants List of descendants made from \code{edge_registry}. Element \code{i} contains all edges descendant from edge \code{i} in the row tree
#' @param colDescendants List of descendants made from \code{edge_registry}. Element \code{j} contains all edges descendant from edge \code{j} in the column tree
#' @param row.tree \code{phylo} class object
#' @examples
incompatible_descendants <- function(rc_ix,rc_tbl,rowDescendants,colDescendants,row.tree){
  row.edg <- rc_tbl[rc_index==rc_ix,row.edge]
  col.edg <- rc_tbl[rc_index==rc_ix,col.edge]
  ### incompatible lineages will have either
  ## the same (or ancestral) col.edges in descendant row.edges
  ## or ancestral row.edges
  
  ## in other words, row_ancs are excluded and only descendant col.edges are allowed in descednant row edges
  
  row_descs <- rowDescendants[[row.edg]]
  col_descs <- colDescendants[[col.edg]]
  row_ancs <- c(edge_ancestors(row.edg,row.tree),row.edg)
  incompatibles <- rc_tbl[row.edge %in% row_ancs,rc_index]
  if (is.null(col_descs)){ ## No possible descendants
    incompatibles <- c(incompatibles,rc_tbl[row.edge %in% row_descs,rc_index])
  } else {
    incompatibles <- c(incompatibles,rc_tbl[row.edge %in% row_descs & !col.edge %in% col_descs,rc_index])
  }
  return(incompatibles)
}