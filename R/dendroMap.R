#' Coarse-grain dataset, X, based on discovered lineage
#' @export
#' @param X dataset input to \code{\link{dendromap}} from which lineages were made
#' @param row.tree \code{phylo} class object
#' @param col.tree \code{phylo} class object
#' @param ncores integer number of cores for parallelization
#' @param maxPval numeric, maximum P-value of \code{row.tree}-x-\code{col.tree} edge pairs to consider
#' @param stepsize granularity of pvals initially input to \code{\link{scan_Fstats}}
#' @examples
dendromap <- function(X,row.tree,col.tree,ncores=NULL,maxPval=0.001,stepsize=10){
  
  # edgeMap, edge_tips for quick descendant calculations -------------------------------------------------------------------
  base::cat('Prepping workspace for edge-based tree traversals')
  rowEdgeMap <- edge_registry(row.tree)
  colEdgeMap <- edge_registry(col.tree)
  row.edges <- 1:Nedge(row.tree)
  col.edges <- 1:Nedge(col.tree)
  rowDescendants <- rowEdgeMap[,list(list(setdiff(seq(edge,edge+n),edge))),by=edge]$V1
  colDescendants <- colEdgeMap[,list(list(setdiff(seq(edge,edge+n),edge))),by=edge]$V1
  colEdgeTips <- edge_tips(col.tree)
  rowEdgeTips <- edge_tips(row.tree)
  
  # edge rc_table -----------------------------------------------------------
  base::cat('\nMaking rc_table')
  rc_table <- make_rc_table(X,row.tree,col.tree,maxPval)
  
  
  # RC_map ------------------------------------------------------------------
  base::cat('\nMaking rc_relations')
  rc_relations <- get_rc_relations(rc_table,rowDescendants,colDescendants)
  
  # simplify rc_table ------------------------------------------------------------------
  base::cat('\nSimplifying rc_table')
  rc_table <- simplify_rc_table(rc_table,rc_relations,rowDescendants,colDescendants)
  rc_relations <- get_rc_relations(rc_table,rowDescendants,colDescendants)
  
  # Scanning F statistics ---------------------------------------------------
  #### set up Pvalues for scanning
  min_pval <-  cbind(rc_table[match(rc_relations$ancestor,rc_index),P],
                     rc_table[match(rc_relations$descendant,rc_index),P]) %>% apply(1,max) %>% min
  rc_relations$max_P <- cbind(rc_table[match(rc_relations$ancestor,rc_index),P],
                              rc_table[match(rc_relations$descendant,rc_index),P]) %>% apply(1,max)
  
  # Prepare cluster -----------------------------------------
  if (!is.null(ncores)){
    base::cat('\nPreparing cluster')
    cl <- parallel::makeCluster(ncores)
    parallel::clusterExport(cl,varlist=c('X','row.tree','col.tree','rc_table','rc_relations',
                                         'rowDescendants','colDescendants','rowEdgeTips','colEdgeTips',
                                         'rowEdgeMap','colEdgeMap'),
                            envir = environment())
    parallel::clusterEvalQ(cl,{library(dendromap)})
  } else {
    cl <- NULL
  }
  
  pvals <- sort(unique(rc_relations$max_P),decreasing = F) ### need an additional trim
  ps <- rc_table[,list(P=min(P),
                       rc_index=rc_index[which.min(P)]),by=row.edge]
  min_P <- rc_relations[ancestor %in% ps$rc_index,max_P] %>% min
  pvals <- pvals[pvals>=min_P]
  scanned_pvals <- pvals[seq(1,length(pvals),by=stepsize)]
  base::cat(paste('Beginning coarse scan of',length(scanned_pvals),'P-value thresholds'))
  
  ### Scanning
  # Fscan <- do.call(scan_Fstats,args = list('pvals'=scanned_pvals,
  #                                          'X'=X,
  #                                          'rc_table'=rc_table,
  #                                          'rc_relations'=rc_relations,
  #                                          'rowDescendants'=rowDescendants,
  #                                          'colDescendants'=colDescendants,
  #                      'rowEdgeTips'=rowEdgeTips,'colEdgeTips'=colEdgeTips,
  #                      'row.tree'=row.tree,
  #                      'cl'=cl),envir = environment())
  Fscan <- scan_Fstats(scanned_pvals,X,rc_table,rc_relations,rowDescendants,colDescendants,
                       rowEdgeTips,colEdgeTips,row.tree,cl=cl)
  # Fscan <- tryCatch(scan_Fstats(scanned_pvals,X,rc_table,rc_relations,rowDescendants,
  #                               rowEdgeTips,colEdgeTips,row.tree,cl=cl),
  #                   error=function(e) e)
  if (!'data.table' %in% class(Fscan)){
    if (!is.null(cl)){
      parallel::stopCluster(cl)
      rm('cl')
      gc()
    }
    stop(Fscan)
  }
  
  
  ### Refining - scan every pval around maximum
  winning_pvals <- Fscan[Fstat==max(Fstat),P_thresh]
  ix <- which(scanned_pvals==min(winning_pvals))
  if (ix<=2){
    refined_pvals <- pvals[pvals<=scanned_pvals[4]]
  } else if (ix>=(length(scanned_pvals)-2)){
    refined_pvals <- pvals[pvals>=scanned_pvals[length(scanned_pvals)-3]]
  } else {
    refined_pvals  <- pvals[pvals>=scanned_pvals[ix-2] & pvals<=scanned_pvals[ix+2]]
  }
  base::cat(paste('Beginning refined scan of',length(refined_pvals),'P-value thresholds'))

  Fscan_refined <- scan_Fstats(refined_pvals,X,rc_table,rc_relations,rowDescendants,
                               colDescendants,rowEdgeTips,colEdgeTips,row.tree,cl=cl)
  # Fscan_refined <- tryCatch(scan_Fstats(refined_pvals,X,rc_table,rc_relations,rowDescendants,
  #                                       rowEdgeTips,colEdgeTips,row.tree,cl=cl),
  #                           error=function(e) e)
  if (!'data.table' %in% class(Fscan_refined)){
    if (!is.null(cl)){
      parallel::stopCluster(cl)
      rm('cl')
      gc()
    }
    stop(Fscan_refined)
  }
  Fscan <- rbind(Fscan,Fscan_refined)
  Fscan <- Fscan[!duplicated(Fscan)]
  
  # Preparing Output ------------------------------------------------------
  
  Lineages <- Fscan[Fstat==max(Fstat),P_thresh] %>% 
    get_lineages_and_stats(rc_table,rc_relations,rowDescendants,
                           colEdgeTips,rowEdgeTips,cl,X,row.tree)
  
  
  if (!is.null(cl)){
    parallel::stopCluster(cl)
    rm('cl')
    gc()
  }
  
  setkey(Lineages,row.edge)
  edge2node <- function(edges,tree) tree$edge[edges,2]
  Lineages[,row.node:=edge2node(row.edge,row.tree)]
  Lineages[,col.node:=edge2node(col.edge,col.tree)]
  Lineages <- Lineages[order(F_stat,decreasing = T)]
  if (!is.null(Lineages$lineage_id)){
    Lineages[,Lineage:=match(lineage_id,unique(lineage_id))]
  }
  setkey(Lineages,F_stat,col.edge)
  object <- list('Lineages'=Lineages,
                 'Data'=X,
                 'row.tree'=row.tree,
                 'col.tree'=col.tree,
                 'colEdgeTips'=colEdgeTips,
                 'rowEdgeTips'=rowEdgeTips,
                 'rowDescendants'=rowDescendants,
                 'colDescendants'=colDescendants,
                 'F_scan'=Fscan)
  class(object) <- 'dendromap'
  return(object)
}