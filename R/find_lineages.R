#' find lineages with RCmap
#' @export
#' @param RCmap see \code{\link{makeRCmap}}
#' @param rc_table see \code{\link{makeRCtable}}
#' @param Row_Descendants named \code{getIndexSets} of all row.nodes
#' @param Col_Descendants named \code{getIndexSets} of all col.nodes
#' @param cl cluster with library \code{dendromap} loaded on each worker
find_lineages <- function(RCmap,rc_table,
                          Row_Descendants,
                          Col_Descendants,cl=NULL){
  nds <- unique(RCmap[terminal==TRUE,descendant])
  base::cat('\nThere are',length(nds),'terminal nodes. Traversing RC tree from all terminal nodes.')
  
  #### RC SEQUENCES = ANT PATHS
  if (is.null(cl)){
    Seqs <- lapply(nds,rc_seqs,RCmap) %>% unlist(recursive=FALSE) %>% unique
  } else {
    Seqs <- parallel::parLapply(cl,nds,rc_seqs,RCmap) %>% unlist(recursive=FALSE) %>% unique
  }
  
  
  #### JOINABILITY
  tbl <- dendromap:::find_joinables(Seqs,rc_table,Row_Descendants,Col_Descendants,cl)

  #### JOINABILITY GRAPH
  seq.indexes <- unique(unlist(tbl[,c('seq1','seq2')]))
  nverts <- length(seq.indexes)
  VertMap <- data.table('seq'=seq.indexes,'vert'=paste('v',seq.indexes,sep='_'),key='vert')
 
  joinables <- tbl[joinability==TRUE]
  nedge <- nrow(joinables)
  joinable.edges <- split(joinables[,c('seq1','seq2')],seq(nrow(joinables))) %>% unlist
  joinable.edges <- VertMap[match(joinable.edges,seq),vert]
  
  base::cat(paste('\nFound',nedge,'joinable RC sequences containing',nverts,' joinable pairs. \n...how hard is it to find maximal cliques in a graph of',length(unique(joinable.edges)),'vertices and',nedge,'edges?'))
  G <- igraph::make_graph(joinable.edges,directed = F)
  
  
  ########## FIND MAXIMAL CLIQUES
  joinable.seqs <- igraph::max_cliques(G,min=2)
  joinable.seqs <- joinable.seqs[order(sapply(joinable.seqs,length),decreasing = F)]
  
  
  lns <- sapply(rev(joinable.seqs),length)
  if (length(joinable.seqs)>10){
    big.lns <- paste(lns[1:5],collapse=',')
    small.lns <- paste(lns[(length(joinable.seqs)-4):length(joinable.seqs)],collapse=',')
  } else {
    small.lns <- paste(lns,collapse=',')
    big.lns <- NULL
  }
  base::cat(paste('\nFound',length(joinable.seqs),'joinable cliques of sizes',big.lns,'...',small.lns))
  base::cat('\nRemoving subsets to find maximal cliques')
  ### below can be parallelized, but now is unncessary b.c. we're using max_cliques?
  # for (i in 1:(length(joinable.seqs)-1)){
  #   check.subset <- any(sapply(joinable.seqs[(i+1):length(joinable.seqs)],
  #                              FUN=function(a,b) all(b %in% a),b=joinable.seqs[[i]]))
  #   if (any(check.subset)){
  #     joinable.seqs[[i]] <- NA
  #   }
  # }
  found.cliques <- !sapply(joinable.seqs,FUN=function(x) all(is.na(x)))
  base::cat(paste('\nThere are',sum(found.cliques),'cliques remaining'))
  joinable.seqs <- joinable.seqs[found.cliques]
  
  ## add disconnected RC's
  all.RC.cliques <- c(as.list(setdiff(VertMap$vert,joinable.edges)),
                     lapply(joinable.seqs,names))
  
  all.RC.cliques <- lapply(all.RC.cliques,FUN=function(v,VertMap) VertMap[match(v,vert),seq],
                        VertMap)
  
  get_rc_index <- function(x,Seqs.=Seqs) unique(unlist(Seqs[x]))
  RCSeqs <- lapply(all.RC.cliques,get_rc_index)
  return(RCSeqs)
}
