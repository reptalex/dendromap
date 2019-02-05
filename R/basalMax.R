#' Find maximum nodes in u along independent nodepaths to root
#' @export
#' @param u vector of length equal to \code{row.tree$Nnode}
#' @param row.tree \code{phylo} class object
#' @param threshold absolute value below which nodes and nodepaths are dropped from consideration. Default is 0
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
#' X <- S$W %*% S$D %*% t(S$V)
#' probs <- binomial(link='logit')$linkinv(X)
#' probs[X==0] <- 0.1
#' P <- rbinom(m*n,1,c(probs))
#' P <- matrix(P,nrow=m,ncol=n,byrow=F)
#' 
#' W <- treeBasis(row.tree)
#' V <- treeBasis(col.tree)
#' U <- t(W) %*% P %*% V
#' 
#' mx <- basalMax(U[,1],row.tree)
#' S$Path[row.node %in% mx$row.node]
basalMax <- function(u,row.tree,threshold=0){
  
  n <- ape::Ntip(row.tree)
  makename <- function(nd) paste('node',nd,sep='_')
  nds <- n+1:row.tree$Nnode
  if (is.null(names(u))){
    names(u) <- makename(nds)
  }
  nms <- names(u)
  
  f <- function(nm) sapply(nm,FUN=function(nm) as.numeric(strsplit(nm,'_')[[1]][2]))
  Maxima <- data.table('row.node'=numeric(round(n/2)),
                       'value'=0)
  k=0
  while (length(u)>0){
    k=k+1
    ix <- which.max(abs(u))
    nd <- f(names(u)[ix])
    
    Maxima$row.node[k] <- nd
    Maxima$value[k] <- u[ix]
    
    rootpath <- ape::nodepath(row.tree,nd,n+1)
    kids <- phangorn::Descendants(row.tree,nd,'all')
    remaining.nodes <- setdiff(f(names(u)),c(rootpath,kids))
    if (length(remaining.nodes)>0){
      u <- u[makename(remaining.nodes)]
    } else {
      u <- NULL
    }
  }
  Maxima <- Maxima[row.node!=0]
  Maxima <- Maxima[abs(value)>=threshold]
  return(Maxima)
}