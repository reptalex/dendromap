#' Find maximum nodes in u along independent nodepaths to root
#' @export
#' @param X Data matrix of number of rows equal to \code{length(row.tree$tip.label)} and number of columns equal to \code{length(col.tree$tip.label)}
#' @param row.tree \code{phylo} class object
#' @param col.tree \code{phylo} class object
#' @param W appropriately named \code{\link{treeBasis}} for \code{row.tree} or subset of that basis. 
#' @param V appropriately named \code{\link{treeBasis}} for \code{col.tree} or subset of that basis. 
#' @param threshold absolute value below which nodes and nodepaths are dropped from consideration. Default is \code{4*sd(log(U^2))}
#' @param min.val minimum value for \eqn{log(U^2)}, where \eqn{U=W'XV}. All values below \code{min.val} will be set to 0
#' @param ncores Number of cores for parallel computation of \code{basalMax}
#' @examples 
#' set.seed(3)
#' m=1e2
#' n=5
#' row.tree <- rtree(m) %>% phytools::force.ultrametric()
#' col.tree <- rtree(n)
#' 
#' S <- treeSim(10,row.tree,col.tree,prob.row=0.8,prob.col=0.8,
#'              col.node = n+1,fix.col.node = T,sd = 1e3,
#'              row.depth.min=2,row.depth.max=3)
#' 
#' X <- makeDataset(S,family='gaussian')
#' Maxima <- basalMaxima(X,row.tree,col.tree,ncores=2)
#' Maxima[,orientation:=sign(value)]
#' 
#' ### code copied from pathViz
#' colmap <- data.table('color'=viridis::viridis(length(unique(Maxima$col.node))),
#'                      'node'=unique(Maxima$col.node))
#' 
#' par(mfrow=c(1,3))
#' plot(row.tree,main='Row tree, raw Gaussian')
#' nodelabels(text=rep(' ',nrow(Maxima)),node=Maxima$row.node,cex=2,
#'            bg = colmap[match(Maxima$col.node,node),color],frame = 'circle')
#' nodelabels(text=rep('.',nrow(S$Paths)),node=S$Paths$row.node,bg='red',frame='circle',cex=.5)
#' 
#' probs <- binomial(link='logit')$linkinv(X)
#' probs[X==0] <- 0.1
#' P <- rbinom(m*n,1,c(probs)) %>% matrix(nrow=m,byrow=F)
#' MaximaP <- basalMaxima(P,row.tree,col.tree,threshold=1.5)
#' plot(row.tree,main='Row tree, Bernoulli data, threshold=1.5')
#' nodelabels(text=rep(' ',nrow(MaximaP)),node=MaximaP$row.node,cex=2,
#'            bg = colmap[match(MaximaP$col.node,node),color],frame = 'circle')
#' nodelabels(text=rep('.',nrow(S$Paths)),node=S$Paths$row.node,bg='red',frame='circle',cex=.5)
#' 
#' MaximaP2 <- basalMaxima(P,row.tree,col.tree,threshold=0.1)
#' plot(row.tree,main='Row tree, Bernoulli data, threshold=0.1')
#' nodelabels(text=rep(' ',nrow(MaximaP2)),node=MaximaP2$row.node,cex=2,
#'            bg = colmap[match(MaximaP2$col.node,node),color],frame = 'circle')
#' nodelabels(text=rep('.',nrow(S$Paths)),node=S$Paths$row.node,bg='red',frame='circle',cex=.5)

basalMaxima <- function(X,row.tree,col.tree=NULL,W=NULL,V=NULL,threshold=NULL,min.val=-10,ncores=NULL){
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
  U <- log(U^2)
  U[U<min.val] <- 0
  if (is.null(threshold)){
    threshold <- mean(U[U>0])+3*sd(U[U>0])
  }

  if (is.null(ncores)){
    Maxima <- apply(U,2,basalMax,row.tree=row.tree,threshold=threshold)
  } else {
    ## Debug
    # getMax <- function(ix,tr,threshold){return(basalMax(U[,ix,drop=F],tr,threshold))}
    # for (ix in colnames(U)){
    #   getMax(ix,row.tree,threshold)
    # }
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
