#' Scan a set of P-value cutoffs, returning F statistics for each
#' 
#' @export
#' @param Pset
#' @param RC_table
#' @param row.nodemap
#' @param col.nodemap
#' @param cl
scanFstats <- function(Pset,RC_table,row.nodemap,col.nodemap,W,V,cl){
  
  Fstats <- numeric(length(Pset))
  for (k in 1:length(Pset)){
    if (k==1){
      rc_table <- RC_table[P<=Pset[k]]
      lineage_table <- findLineageTable(rc_table,row.nodemap,col.nodemap,cl)
    } else {
      ## we only have to recompute lineages from rc_table if they are either descendants/ancestors/same-row.nodes
      ## of row.nodes in rc_table, OR if they are descendants of an ancestor of a new node (sister nodes).
      # We'll use rc_table of new row.nodes - rc_2 - to search for overlaps
      rc_2 <- RC_table[P>Pset[k-1] & P<=Pset[k]] 
      descendants <- lapply(rc_2$row.node,getIndexSets,row.nodemap) %>% unlist %>% unique
      ancestors <- sapply(rc_2$row.node,rootpath,row.tree) %>% unlist %>% c %>% unique
      ancestors_from_prev <- rc_table[row.node %in% ancestors,row.node]
      if (length(ancestors_from_prev)>0){
        cousins <- lapply(ancestors_from_prev,getIndexSets,row.nodemap) %>% unlist %>% unique
      } else {
        cousins <- NULL
      }
      affected_nds <- unique(c(ancestors,descendants,cousins))
      unaffected_lineages <- lineage_table[!row.node %in% affected_nds]
      rc_2 <- rbind(rc_2,rc_table[row.node %in% affected_nds])
      rc_2 <- rc_2[!duplicated(rc_2)]
      new_lineages <- findLineageTable(rc_2,row.nodemap,col.nodemap,cl)
      if (nrow(unaffected_lineages)>0){
        new_lineages$Lineage <- new_lineages$Lineage+max(unaffected_lineages$Lineage)
      }
      lineage_table <- rbind(unaffected_lineages,new_lineages)
      rc_table <- RC_table[P<=Pset[k]]
    }
    Fstats[k] <- getFstat(X,lineage_table,W,V)
  }
  return(Fstats)
}
