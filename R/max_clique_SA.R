#' Simulated annealing max-clique findre
#' @export
#' @param g \code{\link{igraph}} object
#' @param Seqs sequences from \code{\link{rc_seqs}}
#' @param rc_scores scores of \code{rc_index} from \code{\link{makeRCtable}}, with names equal to \code{rc_index} values found in \code{Seqs}
#' @param seq_scores named vector of positive scores. Names must contain all the vertices in \code{g} and match the names of \code{Seqs}
#' @param time.limit time limit, in seconds, for while-loop simultaed annealing
#' @param nreps number of Metropolis-Hastings replicates of randomly drawn cliques to use for initalization
#' @examples 
#' set.seed(1)
#' n=60
#' g <- igraph::random.graph.game(n,0.8)
#' verts <- paste('v_',1:n,sep='')
#' igraph::vertex_attr(g)$name <- verts
#' pscores <- rep(1,length(verts))
#' names(pscores) <- verts
#' 
#' clqs.SA <- max_clique_SA(g,pscores,time.limit=3,nreps=20) 
#' ##time.limit lets us find reasonably large cliques within time.limit
#' 
#' clqs.IG <- igraph::max_cliques(g)  
#' ##this is a good function, but scales poorly with large graphs
#' 
#' lns <- sapply(clqs.IG,length)
#' sort(names(clqs.IG[[which.max(lns)]]))
#' sort(clqs.SA$clique)  
#' ## finds reasonably large cliques, quickly, w/ bonus of weights
#' ## for dendromap, rc-seq cliques will share common root node, so finding quick
#' ## cliques will find the same cophylogenetic lineages as brute-force max_cliques
max_clique_SA <- function(g,Seqs,rc_scores,seq_scores,time.limit=1,nreps=10){
  
  rclique <- function(g,verts,scores){
    if (is.null(verts)){
      verts <- igraph::vertex_attr(g)$name
    }
    if (is.null(scores)){
      scores <- rep(1,length(verts))
      names(scores) <- verts
    }
    A <- igraph::as_adjacency_matrix(g)
    nds <- sample(verts,1,prob=scores)
    remaining <- names(which(A[nds,]>0))
    while (length(remaining)>0){
      nds <- c(nds,sample(remaining,1,prob=scores[match(remaining,verts)]))
      connections <- (matrix(1,ncol=length(nds)) %*% A[nds,])[1,]
      remaining <- names(which(connections==length(nds)))
    }
    return(nds)
  }
  getScore <- function(s,Seqs,rc_scores)    sum(rc_scores[as.character(unique(unlist(Seqs[s])))])
  MH_clique <- function(g,Seqs,nreps,rc_scores,seq_scores,lambda=1){
    verts <- igraph::vertex_attr(g)$name
    nd_list <- lapply(1:nreps,FUN=function(n,g,v,s) rclique(g,v,s),
                      g=g,v=verts,s=seq_scores[verts]^lambda)
    scores <- sapply(nd_list,getScore,Seqs,rc_scores)
    winner <- which.max(scores)
    return(list('clique'=nd_list[[winner]],'score'=scores[[winner]]))
  }
  verts <- igraph::vertex_attr(g)$name
  A <- igraph::as_adjacency_matrix(g)
  one <- matrix(1,ncol=ncol(A))
  
  ### initialize with MH_clique
  clq <- MH_clique(g,Seqs,nreps,rc_scores,seq_scores)
  nds <- clq$clique
  start.time <- Sys.time()
  runtime <- Sys.time()-start.time
  
  i=0
  maxScore <- clq$score
  maxClique <- clq$clique
  while (runtime<=time.limit){
    remaining <- setdiff(verts,nds)
    neighbors <- matrix(1,ncol=length(nds)) %*% A[nds,]
    candidates <- neighbors[,remaining]
    candidates <- candidates[candidates>=(length(nds)-1)]
    if (length(candidates)==0){
      break
    } else if (length(candidates)>1){
      selected.candidate <- sample(candidates,1,prob=seq_scores[names(candidates)])
    } else {
      selected.candidate <- candidates
    }
    candidate.degree <- as.numeric(selected.candidate)
    selected.candidate <- names(selected.candidate)
    selected.candidate.score <- seq_scores[selected.candidate]
    if (candidate.degree<length(nds)){ 
      #some node(s) are not neighbors with our candidate - these are candidates for removal
      removal <- names(which(A[nds,selected.candidate]==0))
      if (length(removal)>1){
        removal <- sample(removal,1,prob=seq_scores[removal])
      }
    } else {
      removal <- sample(nds,1,prob=seq_scores[nds])
    }
    
    temperature <- max(as.numeric(time.limit-runtime),0)
    acceptance.probability <- exp(-(seq_scores[removal]-seq_scores[selected.candidate])/temperature)
    swap <- as.logical(runif(1)<acceptance.probability)
    if (swap){
      nds[nds==removal] <- selected.candidate
    }
    score <- sum(seq_scores[nds])
    if (score>maxScore){
      maxScore <- score
      maxClique <- nds
    }
    runtime <- Sys.time()-start.time
  }
  return(list('clique'=maxClique,'score'=maxScore))
}
