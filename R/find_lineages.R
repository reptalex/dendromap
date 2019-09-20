#' find lineages with RCmap
#' @export
#' @param RCmap see \code{\link{makeRCmap}}
#' @param rc_table see \code{\link{makeRCtable}}
#' @param row.nodemap \code{\link{makeNodeMap}} of row.tree
#' @param col.nodemap \code{\link{makeNodeMap}} of col.tree
find_lineages <- function(RCmap,rc_table,
                          row.nodemap,
                          col.nodemap){
  ### will start with index as descendant and traverse up
  ix <- which(RCmap$terminal)
  base::cat('\nThere are',length(ix),'terminal nodes. Traversing RC tree from all terminal nodes. ')
  
  ### parallelizable
  Seqs <- lapply(ix,rc_seqs,RCmap) %>% unlist(recursive=FALSE) %>% unique
  
  
  n <- length(Seqs)
  tbl <- data.table('seq1'=rep(1:(n-1),times=(n-1):1),key='seq1')
  tbl[,seq2:=(seq1+1):n,by=seq1]
  base::cat(paste('\nFound',length(Seqs),'RC sequences for',nrow(tbl),'pairs. Checking joinability of pairs.'))
  
  
  ### parallelizable, and possibly much more efficient if vectorized w/ data.table
  joinability <- apply(t(tbl),2,FUN=function(x,s,rc,r,c) check_joinable(x[1],x[2],s,rc,r,c),
                       s=Seqs,rc=rc_table,r=row.nodemap,c=col.nodemap)
  tbl[,joinability:=joinability]
  
  
  ### we recursively group joinable sequences into unique index sets...
  ### is joinability transitive? No.
  ### However, it is true that the union of joinable sequences will dominate either sequence
  ### thus we can reduce the number of seqs we check grouping any mutually joinable sets.
  ### these sets will be complete subgraphs in the graph for joinability==TRUE
  ### We can find these with the function igraph::cliques
  joinables <- tbl[joinability==TRUE]
  nedge <- nrow(joinables)
  seq.indexes <- unique(unlist(tbl[,c('seq1','seq2')]))
  nverts <- length(seq.indexes)
  VertMap <- data.table('seq'=seq.indexes,'vert'=paste('v',seq.indexes,sep='_'),key='vert')
  
  base::cat(paste('\nFound',nedge,'joinable RC sequences \ncontaining',nverts,' joinable pairs. \n...how hard is it to find cliques in a graph of',nverts,'vertices and',nedge,'edges?'))
  joinable.edges <- split(joinables[,c('seq1','seq2')],seq(nrow(joinables))) %>% unlist
  joinable.edges <- VertMap[match(joinable.edges,seq),vert]
  
  G <- igraph::make_graph(joinable.edges,directed = F)
  
  base::cat('\nMaking cliques')
  
  ########## FIND MAXIMAL CLIQUES
  ### below can be parallelized by defining subsets (throw out disconnected nodes)
  joinable.seqs <- igraph::cliques(G,min=2) ### now we need to remove elements that are strict subsets of other sets
  joinable.seqs <- joinable.seqs[order(sapply(joinable.seqs,length),decreasing = F)]
  
  
  lns <- sapply(rev(joinable.seqs),length)
  big.lns <- paste(lns[1:5],collapse=',')
  small.lns <- paste(lns[(length(joinable.seqs)-4):length(joinable.seqs)],collapse=',')
  base::cat(paste('\nFound joinable cliques of sizes',big.lns,'...',small.lns))
  base::cat('\nRemoving subsets to find maximal cliques')
  
  ### below can be parallelized
  for (i in 1:(length(joinable.seqs)-1)){
    check.subset <- any(sapply(joinable.seqs[(i+1):length(joinable.seqs)],
                               FUN=function(a,b) all(b %in% a),b=joinable.seqs[[i]]))
    if (any(check.subset)){
      joinable.seqs[[i]] <- NA
    }
  }
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
