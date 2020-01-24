#' find lineages with RCmap
#' @export
#' @param RCmap see \code{\link{makeRCmap}}
#' @param rc_table see \code{\link{makeRCtable}}
#' @param Row_Descendants named \code{getIndexSets} of all row.nodes
#' @param Col_Descendants named \code{getIndexSets} of all col.nodes
#' @param cl cluster with library \code{dendromap} loaded on each worker
#' @param nreps input replicate Metropolis-Hastings simulations to initialize \code{\link{max_clique_SA}}
find_lineages <- function(RCmap,rc_table,
                          Row_Descendants,
                          Col_Descendants,
                          cl,nreps){
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
  
  if (length(Seqs)>1){
    #### JOINABILITY
    tbl <- dendromap:::find_joinables(Seqs,rc_table,Row_Descendants,Col_Descendants,cl)
    if (!is.null(tbl)){
      #### JOINABILITY GRAPH
      seq.indexes <- unique(unlist(tbl[,c('seq1','seq2')]))
      nverts <- length(seq.indexes)
      VertMap <- data.table('seq'=seq.indexes,'vert'=paste('v',seq.indexes,sep='_'),key='vert')
      joinables <- tbl[joinability==TRUE]
      solo_seqs <- setdiff(unique(unlist(tbl[joinability==FALSE,c('seq1','seq2')])),
                             unique(unlist(tbl[joinability==TRUE,c('seq1','seq2')])))
      nedge <- nrow(joinables)
      joinable.edges <- split(joinables[,c('seq1','seq2')],seq(nrow(joinables))) %>% unlist
      joinable.edges <- VertMap[match(joinable.edges,seq),vert]
      
      # base::cat(paste('\nFound',nedge,'joinable RC sequences containing',nverts,' joinable pairs. \n...how hard is it to find maximal cliques in a graph of',length(unique(joinable.edges)),'vertices and',nedge,'edges?'))
      G <- igraph::make_graph(joinable.edges,directed = F)
      
      
      
      ########## FIND MAXIMAL CLIQUES
      sG <- igraph::clusters(G)
      sG.vertices <- lapply(1:sG$no,FUN=function(a,m) names(which(m==a)),m=sG$membership)
      get_rc_index <- function(x,Seqs.=Seqs) unique(unlist(Seqs[x]))
      SubGraphs <- lapply(sG.vertices,igraph::induced_subgraph,graph=G)
      all.RC.cliques <- lapply(SubGraphs,subgraph_to_lineage)
      RCSeqs <- c(lapply(all.RC.cliques,get_rc_index),Seqs[solo_seqs])
    } else {
      RCSeqs <- Seqs
    }
  } else {
    RCSeqs <- Seqs
  }
  return(RCSeqs)
}
