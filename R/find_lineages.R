#' find lineages with RCmap
#' @export
#' @param RCmap see \code{\link{makeRCmap}}
#' @param rc_table see \code{\link{makeRCtable}}
#' @param Row_Descendants named \code{getIndexSets} of all row.nodes
#' @param Col_Descendants named \code{getIndexSets} of all col.nodes
#' @param cl cluster with library \code{dendromap} loaded on each worker
#' @param big_graph_definition integer size of graph (number of vertices) above which \code{\link{max_clique_SA}} will be used to find weighted max clique
#' @param time.limit runtime limit for \code{\link{max_clique_SA}}
#' @param nreps input replicate Metropolis-Hastings simulations to initialize \code{\link{max_clique_SA}}
find_lineages <- function(RCmap,rc_table,
                          Row_Descendants,
                          Col_Descendants,cl,
                          big_graph_definition,
                          time.limit,nreps){
  nds <- unique(RCmap[terminal==TRUE,descendant])
  base::cat('\nThere are',length(nds),'terminal nodes. Traversing RC tree from all terminal nodes.')
  
  seq_score <- function(seq,scores.=scores) sum(scores[as.character(seq)])
  seq_novelty_score <- function(seq,ref,scores.=scores) sum(scores[as.character(setdiff(seq,ref))])
  scaffold_to_rc <- function(scaffold,Seqs.=Seqs) unique(unlist(Seqs[scaffold]))
  
  graph_seqs <- function(g){
    igraph::vertex_attr(g)$name %>% sapply(strsplit,'_') %>%
      sapply('[',2) %>% as.numeric %>% return
  }
  
  most_novel_seq <- function(g,scaffold=NULL,Seqs.=Seqs){
    ix <- graph_seqs(g)
    ix <- setdiff(ix,scaffold)
    ### the graph vertices are indexes of our Seqs. We need to calculate novelty scores for all of those
    ### First, let's get novel rc_indexes in each
    scaffold_rc_indexes <- scaffold_to_rc(scaffold,Seqs)
    novelty_scores <- sapply(Seqs[ix],seq_novelty_score,ref=scaffold_rc_indexes)
    winner <- ix[which.max(novelty_scores)]
    # output <- list(Seqs[[winner]])
    # names(output) <- winner
    return(winner)
  }
  
  remove_incompatible_verts <- function(g,new_winner){
    vert <- paste('v_',new_winner,sep='')
    neighbors <- igraph::neighborhood(g,order = 1,nodes = vert)[[1]]$name
    return(igraph::induced_subgraph(g,neighbors))
  }
  
  check_complete <- function(g,scaffold) all(sort(paste('v_',scaffold,sep='')) ==
                                               sort(igraph::vertex_attr(g)$name))
  
  subgraph_to_lineage <- function(sg){
    scaffold <- NULL
    while(igraph::vcount(sg)>length(scaffold)){
      new_winner <- most_novel_seq(sg,scaffold)
      scaffold <- c(scaffold,new_winner)
      sg <- remove_incompatible_verts(sg,new_winner)
    }
    return(scaffold)
  }
  scores <- -log(rc_table$P)
  names(scores) <- rc_table$rc_index
  
  #### RC SEQUENCES = ANT PATHS
  if (is.null(cl)){
    Seqs <- lapply(nds,rc_seqs,RCmap) %>% unlist(recursive=FALSE) %>% unique
  } else {
    Seqs <- parallel::parLapply(cl,nds,rc_seqs,RCmap) %>% unlist(recursive=FALSE) %>% unique
  }
  names(Seqs) <- paste('v_',1:length(Seqs),sep='')
  
  
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
  
  # base::cat(paste('\nFound',nedge,'joinable RC sequences containing',nverts,' joinable pairs. \n...how hard is it to find maximal cliques in a graph of',length(unique(joinable.edges)),'vertices and',nedge,'edges?'))
  G <- igraph::make_graph(joinable.edges,directed = F)
  
  
  
  ########## FIND MAXIMAL CLIQUES
  sG <- igraph::clusters(G)
  sG.vertices <- lapply(1:sG$no,FUN=function(a,m) names(which(m==a)),m=sG$membership)
  # subgraph_sizes <- sapply(sG.vertices,length)
  # 
  # base::cat(paste('\n--Max subgraph size=',max(subgraph_sizes),
  #                 ' and big_graph_definition=',big_graph_definition,sep=''))
  # if (any(subgraph_sizes>=big_graph_definition)){
  #   SubGraphs <- lapply(sG.vertices,igraph::induced_subgraph,graph=G)
  #   big_graphs <- SubGraphs[subgraph_sizes>=big_graph_definition]
  #   base::cat(paste('\n--Using max_clique_SA on',length(big_graphs),'subgraphs'))
  #   small_graphs <- SubGraphs[subgraph_sizes<big_graph_definition]
  #   rc_scores <- rc_table[,-log(P)]
  #   names(rc_scores) <- rc_table$rc_index
  #   if (is.null(cl)){
  #     seq_scores <- sapply(Seqs,FUN=function(s,rc_table) -sum(log(rc_table[rc_index %in% s,P])),rc_table)
  #   } else {
  #     seq_scores <- parallel::parSapply(cl,Seqs,FUN=function(s,rc_table) -sum(log(rc_table[rc_index %in% s,P])),rc_table)
  #   }
  #   names(seq_scores) <- paste('v_',1:length(Seqs),sep='')
  #   if (is.null(cl)){
  #     big_cliques <- lapply(big_graphs,max_clique_SA,Seqs,rc_scores,seq_scores,time.limit,nreps)
  #   } else {
  #     big_cliques <- parallel::parLapply(cl,big_graphs,max_clique_SA,Seqs,rc_scores,seq_scores,time.limit,nreps)
  #   }
  #   big_cliques <- lapply(big_cliques,getElement,'clique')
  #   small_cliques <- lapply(small_graphs,igraph::max_cliques)
  #   getScore <- function(s,Seqs,rc_scores)    -sum(log(rc_scores[as.character(unique(unlist(Seqs[s])))]))
  #   if (length(small_cliques)>0){
  #     for (i in 1:length(small_cliques)){
  #       if (length(small_cliques[[i]])>1){
  #         scores <- sapply(small_cliques[[i]],getScore,Seqs,rc_scores)
  #         small_cliques[[i]] <- names(small_cliques[[i]][[which.max(scores)]])
  #       } else {
  #         small_cliques[[i]] <- names(small_cliques[[i]][[1]])
  #       }
  #     }
  #   }
  #   joinable.seqs <- c(big_cliques,small_cliques)
  # } else {
  #   joinable.seqs <- igraph::max_cliques(G,min=2)
  # }
  # joinable.seqs <- joinable.seqs[order(sapply(joinable.seqs,length),decreasing = F)]
  # 
  # 
  # lns <- sapply(rev(joinable.seqs),length)
  # if (length(joinable.seqs)>10){
  #   big.lns <- paste(lns[1:5],collapse=',')
  #   small.lns <- paste(lns[(length(joinable.seqs)-4):length(joinable.seqs)],collapse=',')
  # } else {
  #   small.lns <- paste(lns,collapse=',')
  #   big.lns <- NULL
  # }
  # base::cat(paste('\nFound',length(joinable.seqs),'joinable cliques of sizes',big.lns,'...',small.lns))
  # # base::cat('\nRemoving subsets to find maximal cliques')
  # ### below can be parallelized, but now is unncessary b.c. we're using max_cliques?
  # # for (i in 1:(length(joinable.seqs)-1)){
  # #   check.subset <- any(sapply(joinable.seqs[(i+1):length(joinable.seqs)],
  # #                              FUN=function(a,b) all(b %in% a),b=joinable.seqs[[i]]))
  # #   if (any(check.subset)){
  # #     joinable.seqs[[i]] <- NA
  # #   }
  # # }
  # found.cliques <- !sapply(joinable.seqs,FUN=function(x) all(is.na(x)))
  # joinable.seqs <- joinable.seqs[found.cliques]
  # 
  # ## add disconnected RC's
  # all.RC.cliques <- c(as.list(setdiff(VertMap$vert,joinable.edges)),
  #                    lapply(joinable.seqs,names))
  # 
  # all.RC.cliques <- lapply(all.RC.cliques,FUN=function(v,VertMap) VertMap[match(v,vert),seq],
  #                       VertMap)
  # 
  get_rc_index <- function(x,Seqs.=Seqs) unique(unlist(Seqs[x]))
  SubGraphs <- lapply(sG.vertices,igraph::induced_subgraph,graph=G)
  all.RC.cliques <- lapply(SubGraphs,subgraph_to_lineage)
  RCSeqs <- lapply(all.RC.cliques,get_rc_index)
  return(RCSeqs)
}
