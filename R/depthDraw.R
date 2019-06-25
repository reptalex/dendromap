#' draw node pairs for treeMap
#' 
#' @param d1 depths on row tree
#' @param d2 depths on col tree
#' @param lambda propensities of column tree nodes
#' @param decay.rate exponential rate of decay for past nodes. Probability of draw will be exp(-decay.rate*dT) where dT is the time from column node to row node.
depthDraw <- function(d1,d2,lambda,decay.rate=0){
  
  D <- data.table(expand.grid('row'=1:length(d1),'col.node'=1:length(d2)))
  D[,lambda:=lambda[col.node]]
  D[,dT:=d1[row]-d2[col.node]]
  D[dT<0,dT:=0]
  
  sampler <- function(col.node,dT,lambda,decay.rate){
    if (all(lambda==-Inf)){
      return(as.integer(0))
    } else {
      happens=c(T,F)[rbinom(1,1,prob = exp(-sum(lambda[!lambda==-Inf])))+1] & !all(dT==0)
      if (happens){
        col.node <- col.node[dT>0]
        lambda <- lambda[dT>0]
        dT <- dT[dT>0]
        if (any(is.infinite(lambda))){
          
          if (all(lambda[is.infinite(lambda)]<0)){
            vals <- !is.infinite(lambda)
            lambda <- lambda[vals]
            if (length(lambda)==0){
              cs <- 0
            } else {
              cs <- sample(col.node[ix],1,prob=lambda*exp(-decay.rate*dT))
            }
          } else {  ## some positive infinity vals
            vals <- is.infinite(lambda) & lambda>0
            if (sum(vals)==1){
              cs <- col.node[vals]
            } else {
              cs <- sample(col.node[vals],1,prob=rep(1,sum(vals))*exp(-decay.rate*dT[vals]))
            }
          }
        } else {
          cs <- sample(col.node,1,prob=lambda*exp(-decay.rate*dT))
        }
        return(as.integer(cs))
      } else {
        return(as.integer(0))
      }
    }
  }
    
  D[,list(col=sampler(col.node,dT,lambda,decay.rate)),by=row] %>%
  return()
}
