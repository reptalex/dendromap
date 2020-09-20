#' Remove all (r,c) such that P(r,c)>P(r',c) for some descendant r' of r 
#' @export
#' @param rc_tbl see \code{\link{make_rc_table}}
#' @param rowDescendants List of descendants made from \code{edge_registry}. Element \code{i} contains all edges descendant from edge \code{i} in the row tree
#' @examples
basal_dominance_filter <- function(rc_tbl,rowDescendants){
  ### Filtering ### 
  ## We'd like to do this as few times as possible.
  ## One option is to start with basal indexes of rc_table for the min P-value cutoff
  ## and then use that filtered rc_table's basal-index shards for assemble_lineage input
  ## When adding more rc-indexes with higher P-values, 
  ## we can simply check whether/not they, too, need to be excluded.
  
  ## RULES OF FILTERING:
  ## remove col.edges with more significant values among descendants (implies signal lost with inclusion of sister taxa)
  ### Implications:
  ## addition of P-values with higher values will not upset a lower P-value's filtering. 
  setkey(rc_tbl,row.edge,P)
  row.edges <- rc_tbl[,unique(row.edge)]
  for (ee in row.edges){
    dum <- rc_tbl[row.edge==ee,c('col.edge','P','rc_index')]
    if (nrow(dum)>0 & length(rowDescendants[[ee]])>0){
      descs_tbl <- rc_tbl[row.edge %in% rowDescendants[[ee]] & col.edge %in% dum$col.edge]
      descs <- descs_tbl[,list(Pmin=min(P)),by=col.edge]
      setkey(dum,col.edge)
      setkey(descs,col.edge) ## this may not contain all col.edges in dum
      dd <- descs[dum]
      removers <- dd[P>Pmin,rc_index]
      rc_tbl <- rc_tbl[!rc_index %in% removers]
    }
  }
  return(rc_tbl)
}