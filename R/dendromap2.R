
#' dendromap - find cophylogenetic patterns in a dataset
#' @export
#' @param X matrix whose rownames are in \code{row.tree$tip.label} and colnames are in \code{col.tree$tip.label}
#' @param row.tree phylo class object. Polytomies will be ignored.
#' @param col.tree phylo class object. Polytomies will be ignored.
#' @param ncores If not NULL, then integer specifying number of cores to use for parallelizable steps
#' @param max_Pval maximum P-value for consideration - will scan all P<max_Pval for optimal lineage.
#' @param nP maximum number of P-values between the minimum observed P-value and the \code{max_Pval} specified to scan.
#' @param W optional \code{treeBasis(row.tree)} - must have colnames in the format of e.g. "node_51" for node 51.
#' @param V optional \code{treeBasis(col.tree)} - must have colnames in the format of e.g. "node_51" for node 51.
#' @param n_sim optional number of null datasets to simulate (via row & column shuffling) in order to generate approximate P-values in \code{\link{makeRCtable}}
#' @examples
#' library(dendromap)
#' set.seed(3)
#' m=1e3
#' n=30
#' row.tree <- rtree(m) %>% phytools::force.ultrametric()
#' col.tree <- rtree(n)
#' S <- treeSim(5,row.tree,col.tree,row.depth.min=2,row.depth.max=3,col.node=n+1,fix.col.node=T) 
#' eta <- S$W %*% (10*sign(S$D)) %*% t(S$V)
#' X <- eta+matrix(rnorm(m*n),nrow=m)
#' clrinv <- function(x) exp(x)/sum(exp(x))
#' rmlt <- function(p,lambda=5e3) rmultinom(1,rpois(1,lambda),prob = p)
#' N <- apply(X,2,clrinv) %>% apply(2,rmlt)
#' rownames(N) <- row.tree$tip.label
#' colnames(N) <- col.tree$tip.label
#' dm <- dendromap(N,row.tree,col.tree,Pval_threshold=0.05)
#' ## can use multiple cores for parallelization. Will speed-up large datasets.
#' ## big graphs have to be handled with an alternative max_clique algorithm: max_clique_SA
#' ## dm2 <- dendromap(N,row.tree,col.tree,W=S$W,V=S$V,ncores=2,Pval_threshold=0.2)
#' #Since they've already been computed, inputting the matrices W, V saves time.
#' 
#' dendromap:::print.dendromap(S)
#' dendromap:::print.dendromap(dm)
#' 
#' dendromap:::plot.dendromap(S)
#' dendromap:::plot.dendromap(dm)
#' 
#' S$Lineages[,rc:=paste(row.node,col.node,sep='_')]
#' dm$Lineages[,rc:=paste(row.node,col.node,sep='_')]
#' sum(S$Lineages$rc %in% dm$Lineages$rc)/nrow(S$Lineages)      ### Probability of a true rc being ID'd
#' sum(dm$Lineages$rc %in% S$Lineages$rc)/nrow(dm$Lineages)     ### Probability of ID'd rc being true positive

dendromap2 <- function(X,row.tree,col.tree,ncores=NULL,
                      max_Pval=0.01,nP=100,W=NULL,V=NULL,n_sim=NULL,
                      estimate_runtime=FALSE){
  
  base::cat(paste('Checking Data and tree compatibility'))
  ### Align dataset to trees
  if (!all(rownames(X) %in% row.tree$tip.label)){
    stop('There are rownames(X) not in row.tree$tip.label')
  } else {
    if (!all(row.tree$tip.label %in% rownames(X))){
      row.tree <- ape::drop.tip(row.tree,setdiff(rownames(X),row.tree$tip.label))
    }
    X <- X[row.tree$tip.label,]
  }
  if (!all(colnames(X) %in% col.tree$tip.label)){
    stop('There are colnames(X) not in col.tree$tip.label')
  } else {
    if (!all(col.tree$tip.label %in% colnames(X))){
      col.tree <- ape::drop.tip(col.tree,setdiff(colnames(X),col.tree$tip.label))
    }
    X <- X[,col.tree$tip.label]
  }
  if (!is.null(ncores)){
    cl <- parallel::makeCluster(ncores)
    parallel::clusterEvalQ(cl,library(dendromap))
  } else {
    cl <- NULL
  }
  if (is.null(W)){
    base::cat(paste('Making treeBasis for row.tree'))
    W <- treeBasis(row.tree)
  }
  if (is.null(V)){
    base::cat(paste('Making treeBasis for col.tree'))
    V <- treeBasis(col.tree)
  }
  
  base::cat(paste('\nMaking Nodemaps'))
  row.nodemap <- dendromap:::makeNodeMap(row.tree)
  col.nodemap <- dendromap:::makeNodeMap(col.tree)
  
  base::cat(paste('\nMaking RC table with',n_sim,'null simulations'))
  RC_table <- makeRCtable(X,row.tree,col.tree,W,V,n_sim)
  Pset <- unique(RC_table[P<=max_Pval,P]) %>% sort(decreasing=F)
  if (length(Pset)>nP){
    Pset <- seq(min(Pset),max_Pval,length.out=nP)
  }
  
  if (estimate_runtime){
    t1 <- Sys.time()
    findLineageTable(RC_table[P<=max(Pset)],row.nodemap,col.nodemap,cl)
    t2 <- Sys.time()
    base::cat(paste('\nCompute time for the largest set of nodes took',signif(t2-t1,3),'seconds. \nEstimated upper-bound time of completion:',Sys.time()+length(Pset)*(t2-t1)))
  }
  
  # return(RC_table)
  Fstats <- scanFstats(Pset,RC_table,row.nodemap,col.nodemap,W,V,cl)
  k <- which.max(Fstats)
  OptimalP <- Pset[k]
  output <- findLineageTable(RC_table[P<=OptimalP],row.nodemap,col.nodemap,cl)
  
  if (!is.null(cl)){
    parallel::stopCluster(cl)
    rm('cl')
    gc()
  }
  output <- list('Lineages'=output,
                 'Data'=X,
                 'row.tree'=row.tree,
                 'col.tree'=col.tree,
                 'Pvalues_scanned'=Pset,
                 'Fstats'=Fstats,
                 'winning_Pval_index'=k)
  class(output) <- 'dendromap'
  return(output)
}
