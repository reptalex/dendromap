#' find lineages with RCmap
#' @param RCmap see \code{\link{makeRCmap}}
#' @param rc_table see \code{\link{makeRCtable}}
#' @param row.nodemap \code{\link{makeNodeMap}} of row.tree
#' @param col.nodemap \code{\link{makeNodeMap}} of col.tree
find_lineages <- function(RCmap,rc_table,
                         row.nodemap.=row.nodemap,
                         col.nodemap.=col.nodemap){
  ### will start with index as descendant and traverse up
  ix <- which(RCmap$terminal)
  Seqs <- lapply(ix,rc_seqs,RCmap) %>% unlist(recursive=FALSE) %>% unique
  n <- length(Seqs)
  
  
  tbl <- data.table('seq1'=rep(1:(n-1),times=(n-1):1),key='seq1')
  tbl[,seq2:=(seq1+1):n,by=seq1]
  
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
  joinable.edges <- split(joinables[,c(1,2)],seq(nrow(joinables))) %>% unlist
  G <- igraph::make_graph(joinable.edges,directed = F)
  joinable.seqs <- igraph::cliques(G) ### now we need to remove elements that are strict subsets of other sets
  joinable.seqs <- joinable.seqs[order(sapply(joinable.seqs,length),decreasing = F)]
  for (i in 1:(length(joinable.seqs)-1)){
    check.subset <- any(sapply(joinable.seqs[(i+1):length(joinable.seqs)],
                               FUN=function(a,b) all(b %in% a),b=joinable.seqs[[i]]))
    if (any(check.subset)){
      joinable.seqs[[i]] <- NA
    }
  }
  found.cliques <- !sapply(joinable.seqs,FUN=function(x) all(is.na(x)))
  joinable.seqs <- joinable.seqs[found.cliques]
  joinable.seqs <- c(setdiff(unique(unlist(tbl[,c(1,2)])),joinable.seqs),
                     joinable.seqs)
  get_rc_index <- function(x,Seqs.=Seqs) unique(unlist(Seqs[x]))
  Seqs <- lapply(joinable.seqs,get_rc_index)
  return(Seqs)
}