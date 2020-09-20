scanFstats <- function(pvals,rc_table.=rc_table,RCmap.=RCmap,cl=NULL){
  ix_prev <- NULL
  Fstats <- rep(NA,length(pvals))
  if (!is.null(cl)){
    ### check cluster?
  }
  
  for (i in 1:length(pvals)){
    p_thresh=pvals[i]
    ix_thresh=rc_table[P<=p_thresh,rc_index]
    new_ix <- setdiff(ix_thresh,ix_prev)
    rct <- rc_table[P<=p_thresh]
    rcm <- RCmap[max_P<=p_thresh]
    
    if (i==1){
      basal_indexes <- rcm[ancestor %in% ix_thresh & descendant %in% ix_thresh,setdiff(unique(ancestor),unique(descendant))]
      desc_count <- rcm[ancestor %in% ix_thresh & descendant %in% ix_thresh,list(n=.N),by=ancestor]
      basal_ixs_with_descendants <- intersect(basal_indexes,desc_count[n>1,ancestor])
      if (is.null(cl)){
        Lineages <- lapply(basal_ixs_with_descendants,getLins,p_thresh=p_thresh) %>% rbindlist
      } else {
        Lineages <- parallel::parLapply(cl,basal_ixs_with_descendants,getLins,p_thresh=p_thresh) %>% rbindlist
      }
      Lineages[,lineage_size:=.N,by=lineage_id]
      Lineages <- Lineages[lineage_size>1]
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
      new_descs <- Row_Descendants[new_row.edges]
      new_anc_edges <- sapply(new_descs,FUN=function(a,b) any (b %in% a), b=prev_row.edges) %>% new_row.edges[.]
      
      ## find previous edges ancestral to new edges
      prev_descs <- Row_Descendants[prev_row.edges]
      prev_anc_edges <- sapply(prev_descs,FUN=function(a,b) any (b %in% a), b=new_row.edges) %>% prev_row.edges[.]
      
      ## need to recompute all edges which are either equal, ancestral, or descendant to new edges
      row_path_edges <- c(intersect(prev_row.edges,new_row.edges),
                          prev_anc_edges,
                          new_anc_edges) %>% unique
      
      row_path_ix <- rct[row.edge %in% row_path_edges,unique(rc_index)]
      ## NOTE: Since clean_sisters can occasionally introduce a new row.edge not found
      ##       we'll remove all lineage_id's containing row_edges in our row_path_edges set above
      Old_Lineages <- Lineages[!row.edge %in% row_path_edges]
      
      basal_ix <- rcm[ancestor %in% row_path_ix & descendant %in% ix_thresh,setdiff(unique(ancestor),unique(descendant))]
      desc_count <- rcm[ancestor %in% row_path_ix & descendant %in% ix_thresh,list(n=.N),by=ancestor]
      basal_ixs_with_descendants <- intersect(basal_ix,desc_count[n>1,ancestor])
      if (is.null(cl)){
        Lineages <- lapply(basal_ixs_with_descendants,getLins,p_thresh=p_thresh) %>% rbindlist
      } else {
        Lineages <- parallel::parLapply(cl,basal_ixs_with_descendants,getLins,p_thresh=p_thresh) %>% rbindlist
      }
      Lineages[,lineage_size:=.N,by=lineage_id]
      Lineages <- Lineages[lineage_size>1]
      Lineages <- rbind(Old_Lineages,Lineages)
    }
    
    
    #### need to parallelize lineage_stats
    if (nrow(Lineages)>0){
      stats <- lineage_stats(Lineages,X,colEdgeTips,rowEdgeTips,cl=cl)
      stats <- stats[order(F_stat,decreasing = T)]
      stats <-  filter_stats(stats,Lineages,rct,rcm,Row_Descendants)
      Lineages <- Lineages[lineage_id %in% stats$lineage_id]
    }
    if (nrow(Lineages)>0){
      Fstats[i] <- global_Fstat(Lineages,X,colEdgeTips,rowEdgeTips)
    }
    
    ix_prev <- ix_thresh  ## we only have to recompute new_ix that affect our Lineages table
  }
  
  return(data.table('P_thresh'=pvals,'Fstat'=Fstats))
}
