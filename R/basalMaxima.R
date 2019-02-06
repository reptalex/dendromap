#' Find maximum nodes in u along independent nodepaths to root
#' @export
#' @param X Data matrix of number of rows equal to \code{length(row.tree$tip.label)} and number of columns equal to \code{length(col.tree$tip.label)}
#' @param row.tree \code{phylo} class object
#' @param col.tree \code{phylo} class object
#' @param W appropriately named \code{\link{treeBasis}} for \code{row.tree} or subset of that basis. 
#' @param V appropriately named \code{\link{treeBasis}} for \code{col.tree} or subset of that basis. 
#' @param threshold absolute value below which nodes and nodepaths are dropped from consideration. Default is 0
#' @param ncores Number of cores for parallel computation of \code{basalMax}
#' @examples 
#' set.seed(3)
#' m=1e3
#' n=30
#' row.tree <- rtree(m) %>% phytools::force.ultrametric()
#' col.tree <- rtree(n)
#' 
#' S <- treeSim(10,row.tree,col.tree,prob.row=0.7,prob.col=0.8,
#'              col.node = n+1,fix.col.node = T,sd = 1e3,
#'              row.depth.min=2,row.depth.max=3)
#' 
#' X <- makeDataset(S,family='gaussian')
#' Maxima <- basalMaxima(X,row.tree,col.tree,ncores=2)
#' Maxima[,orientation:=sign(value)]
#' 
#' ### code copied from pathViz
#' colmap <- viridis::viridis(nrow(Maxima))
#' par(mfrow=c(1,1))
#' plot(row.tree,main='Row tree')
#' nodelabels(text=rep(' ',nrow(Maxima)),node=Maxima$row.node,bg = colmap,frame = 'circle')
#' nodelabels(text=rep('.',nrow(S$Paths)),node=S$Paths$row.node,bg='red',frame='circle',cex=.5)

basalMaxima <- function(X,row.tree,col.tree=NULL,W=NULL,V=NULL,threshold=5,ncores=NULL){
  if (is.null(col.tree) & is.null(V)){
    stop('Must input either col.tree or V')
  }
  if (is.null(W)){
    W <- treeBasis(row.tree)
  }
  if (is.null(V)){
    V <- treeBasis(col.tree)
  }
  
  U <- t(W) %*% X %*% V
  if (is.null(ncores)){
    Maxima <- apply(U,2,basalMax,row.tree=S$row.tree,threshold=threshold)
  } else {
    cl <- parallel::makeCluster(ncores)
    getMax <- function(ix,tr,threshold){return(basalMax(U[,ix,drop=F],tr,threshold))}
    parallel::clusterExport(cl,varlist=c('basalMax','U','getMax'),envir = environment())
    parallel::clusterEvalQ(cl,expr=library(data.table))
    Maxima <- NULL
    Maxima <- tryCatch(parallel::parLapply(cl,colnames(U),getMax,row.tree,threshold),
                       error=function(e) as.character(e))
    parallel::stopCluster(cl)
    if (class(Maxima)=='character'){
      stop(paste('error in parLapply(cl,colnames(U),getMax,row.tree,threshold):',Maxima))
    } 
  }
  for (i in 1:length(Maxima)){
    if (nrow(Maxima[[i]])>0){
      nd <- as.numeric(strsplit(colnames(U)[i],'_')[[1]][2])
      Maxima[[i]][,col.node:=nd]
    }
  }
  nr <- sapply(Maxima,nrow)
  Maxima <- Maxima[nr>0]
  Maxima <- data.table::rbindlist(Maxima)
  return(Maxima)
}
