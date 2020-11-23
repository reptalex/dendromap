#' Get lineages with all lineage F-statistics
#' @export
#' @param p_thresh Positive numeric, threshold P-value for inclusion of rc-pairs
#' @param rc_table see \code{\link{make_rc_table}}
#' @param rc_relations see \code{\link{get_rc_relations}}
#' @param rowDescendants List of descendants made from \code{edge_registry}. Element \code{i} contains all edges descendant from edge \code{i} in the row tree
#' @param colDescendants List of descendants made from \code{edge_registry}. Element \code{j} contains all edges descendant from edge \code{j} in the column tree
#' @param colEdgeTips \code{\link{edge_tips}} for col.tree
#' @param rowEdgeTips \code{\link{edge_tips}} for row.tree
#' @param cl cluster initialized internally with \code{\link{dendromap}} 
#' @param X matrix. Dataset from which F-statistics are computed
#' @param row.tree phylo class object, with tips corresponding to rows in \code{X}
#' @examples
get_lineages_and_stats <- function(p_thresh,rc_table,rc_relations,rowDescendants,
                                   colEdgeTips,rowEdgeTips,cl=NULL,X,row.tree){
  ix_thresh=rc_table[P<=p_thresh,rc_index]
  rct <- rc_table[rc_index %in% ix_thresh]
  rcm <- rc_relations[max_P<=p_thresh]
  basal_indexes <- rcm[ancestor %in% ix_thresh & descendant %in% ix_thresh,setdiff(unique(ancestor),unique(descendant))]
  desc_count <- rcm[ancestor %in% ix_thresh & descendant %in% ix_thresh,list(n=.N),by=ancestor]
  basal_ixs_with_descendants <- intersect(basal_indexes,desc_count[n>1,ancestor])
  if (is.null(cl)){
    Lineages <- lapply(basal_ixs_with_descendants,get_lineage,p_thresh,
                       rc_relations=rc_relations,rc_table=rc_table,row.tree=row.tree,
                       rowDescendants=rowDescendants,colDescendants=colDescendants) %>% rbindlist
  } else {
    Lineages <- parallel::parLapply(cl,basal_ixs_with_descendants,get_lineage_par,p_thresh) %>% rbindlist
  }
  Lineages[,lineage_size:=.N,by=lineage_id]
  Lineages <- Lineages[lineage_size>1]
  stats <- lineage_stats(Lineages,X,colEdgeTips,rowEdgeTips)
  stats <- stats[order(F_stat,decreasing = T)]
  stats <-  filter_stats(stats,Lineages,rct,rcm,rowDescendants,row.tree)
  Lineages <- Lineages[lineage_id %in% stats$lineage_id]
  
  setkey(Lineages,lineage_id,row.edge)
  setkey(stats,lineage_id)
  
  Lineages <- stats[,c('lineage_id','F_stat')][Lineages]
  return(Lineages)
}
