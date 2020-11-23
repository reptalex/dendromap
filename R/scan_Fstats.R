#' Scan P-value thresholds, computing global F-statistic for lineages obtained at each step
#' @export
#' @param pvals numeric vector containing values between 0 and 1, sorted from smallest to largest
#' @param X dataset input to \code{\link{dendromap}}
#' @param rc_table see \code{\link{make_rc_table}}
#' @param rc_relations see \code{\link{get_rc_relations}}
#' @param rowDescendants List of descendants made from \code{edge_registry}. Element \code{i} contains all edges descendant from edge \code{i} in the row tree
#' @param colDescendants List of descendants made from \code{edge_registry}. Element \code{i} contains all edges descendant from edge \code{i} in the col tree
#' @param rowEdgeTips \code{\link{edge_tips}} for row.tree
#' @param colEdgeTips \code{\link{edge_tips}} for col.tree
#' @param row.tree phylo class object
#' @param cl cluster initialized internally with \code{\link{dendromap}} 
#' @examples
scan_Fstats = function(pvals,X,rc_table,rc_relations,
                       rowDescendants,colDescendants,rowEdgeTips,colEdgeTips,row.tree,cl=NULL){
  ix_prev <- NULL
  Fstats <- rep(NA,length(pvals))
  start_time <- Sys.time()
  Lineages <- data.table()
  
  # if (is.null(cl)){
  #   dm_env <- new.env()
  #   dm_env$rc_table <- rc_table
  #   dm_env$rc_relations <- rc_relations
  #   dm_env$rowDescendants <- rowDescendants
  #   dm_env$colDescendants <- colDescendants
  #   dm_env$rowEdgeTips <- rowEdgeTips
  #   dm_env$colEdgeTips <- colEdgeTips
  #   dm_env$row.tree <- row.tree
  # }
  
  for (i in 1:length(pvals)){
    p_thresh=pvals[i]
    ix_thresh=rc_table[P<=p_thresh,rc_index]
    new_ix <- setdiff(ix_thresh,ix_prev)
    rct <- rc_table[P<=p_thresh]
    rcm <- rc_relations[max_P<=p_thresh]
    
    if (nrow(Lineages)==0){
      basal_indexes <- rcm[ancestor %in% ix_thresh & 
                             descendant %in% ix_thresh,
                           setdiff(unique(ancestor),unique(descendant))]
      desc_count <- rcm[ancestor %in% ix_thresh & descendant %in% ix_thresh,
                        list(n=.N),by=ancestor]
      basal_ixs_with_descendants <- intersect(basal_indexes,desc_count[n>=1,ancestor])
      if (length(basal_ixs_with_descendants)>0){
        if (is.null(cl)){
          # Lineages <- with(environment(),
          #                   lapply(basal_ixs_with_descendants,get_lineage,p_thresh=p_thresh,
          #                   rc_relations=rc_relations,rc_table=rc_table,row.tree=row.tree,
          #                   rowDescendants=rowDescendants,colDescendants=colDescendants))
          # Lineages <- lapply(basal_ixs_with_descendants,get_lineage,p_thresh=p_thresh,
          #                    rc_relations=rc_relations,rc_table=rc_table,row.tree=row.tree,
          #                    rowDescendants=rowDescendants,colDescendants=colDescendants) %>%
          Lineages <- lapply(basal_ixs_with_descendants,get_lineage,p_thresh=p_thresh,rc_table=rc_table) %>%
            rbindlist
        } else {
          Lineages <- parallel::parLapply(cl,basal_ixs_with_descendants,get_lineage,p_thresh=p_thresh) %>% 
            rbindlist
        }
        
        Lineages[,lineage_size:=.N,by=lineage_id]
        Lineages <- Lineages[lineage_size>1]
      }
    } else {
      ### here, we note that our simplify_rc_table function ensured that any new ix will:
      ### if a descendant row edge from a basal_ix, it cannot have same col.edge as prior basal_ix 
      ##   - this is because P(r,c) < P(r',c) for r' descendant of r <--> basal dominance.
      
      ### Which basal_indexes to we (re)compute?
      ### - ancestors/descendants of new_ix
      ### - note: we may even produce conflicting lineages, e.g. Rheas Gondwanan and Rheas Australian
      ### we DON'T compute: basal_ix NOT on row-tree path (root-tips) from new_ix.
      ### so... 
      ### (1) row_path_ix: Set of ix_thresh on row.tree root-tip path of new_ix
      ### (2) remove all Lineages containing row_path_ix
      ### (3) find basal_ix among row_path_ix
      ### (4) compute Lineages for basal_ix
      ### (5) rbind with old lineages
      ### (6) stats-->filter-->global F stat
      
      prev_row.edges <- rct[rc_index %in% ix_prev,unique(row.edge)]
      new_row.edges <- rct[rc_index %in% new_ix,unique(row.edge)]
      
      ## find new adges ancestral to previous edges
      new_descs <- rowDescendants[new_row.edges]
      new_anc_edges <- sapply(new_descs,FUN=function(a,b) any (b %in% a), b=prev_row.edges) %>%
        new_row.edges[.]
      
      ## find previous edges ancestral to new edges
      prev_descs <- rowDescendants[prev_row.edges]
      prev_anc_edges <- sapply(prev_descs,FUN=function(a,b) any (b %in% a), b=new_row.edges) %>% 
        prev_row.edges[.]
      
      ## need to recompute all edges which are either equal, ancestral, or descendant to new edges
      row_path_edges <- c(intersect(prev_row.edges,new_row.edges),
                          prev_anc_edges,
                          new_anc_edges) %>% unique
      
      row_path_ix <- rct[row.edge %in% row_path_edges,unique(rc_index)]
      ## NOTE: Since merge_twin_sisters can occasionally introduce a new row.edge not found
      ##       we'll remove all lineage_id's containing row_edges in our row_path_edges set above
      Old_Lineages <- Lineages[!row.edge %in% row_path_edges]
      
      basal_ix <- rcm[ancestor %in% row_path_ix & descendant %in% ix_thresh,
                      setdiff(unique(ancestor),unique(descendant))]
      desc_count <- rcm[ancestor %in% row_path_ix & descendant %in% ix_thresh,
                        list(n=.N),by=ancestor]
      basal_ixs_with_descendants <- intersect(basal_ix,desc_count[n>=1,ancestor])
      if (length(basal_ixs_with_descendants)>0){
        if (is.null(cl)){
          # Lineages <- lapply(basal_ixs_with_descendants,get_lineage,p_thresh=p_thresh,
          #                    rc_relations=rc_relations,rc_table=rc_table,row.tree=row.tree,
          #                    rowDescendants=rowDescendants,colDescendants=colDescendants) %>% 
          #   rbindlist
          Lineages <- lapply(basal_ixs_with_descendants,get_lineage,p_thresh=p_thresh) %>%
            rbindlist
        } else {
          Lineages <- parallel::parLapply(cl,basal_ixs_with_descendants,get_lineage,p_thresh=p_thresh) %>% 
            rbindlist
        }
        
        Lineages[,lineage_size:=.N,by=lineage_id]
        Lineages <- Lineages[lineage_size>1]
        Lineages <- rbind(Old_Lineages,Lineages)
      }
    }
    
    
    if (nrow(Lineages)>0){
      stats <- lineage_stats(Lineages,X,colEdgeTips,rowEdgeTips,cl=cl)
      stats <- stats[order(F_stat,decreasing = T)]
      stats <-  filter_stats(stats,Lineages,rct,rcm,rowDescendants,row.tree)
      Lineages <- Lineages[lineage_id %in% stats$lineage_id]
    }
    if (nrow(Lineages)>0){
      Fstats[i] <- global_Fstat(Lineages,X,colEdgeTips,rowEdgeTips)
    }
    
    ix_prev <- ix_thresh  ## we only have to recompute new_ix that affect our Lineages table
    tm2 <- Sys.time()
    time.elapsed <- signif(difftime(tm2,start_time,units = 'mins'),3)
    GUI.notification <- paste('\r',i,'P thresholds out of',length(pvals),
                              'scanned in',time.elapsed,'minutes.    ')
    GUI.notification <- paste(GUI.notification,'Estimated time of completion for this step:',
                              as.character(start_time+difftime(tm2,start_time)*length(pvals)/i),
                              '  \r')
    base::cat(GUI.notification)
    utils::flush.console()
  }
  return(data.table('P_thresh'=pvals,'Fstat'=Fstats))
}
