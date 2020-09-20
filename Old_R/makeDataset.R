#' make dataset from treemap object
#' @export
#' @param Path output from \code{\link{treeMap}}
#' @param row.tree \code{phylo} class object used to make \code{Path}
#' @param col.tree \code{phylo} class object used to make \code{Path}
#' @param d numeric vector of length equal to \code{nrow(Path)}
#' @param sd standard deviation for simulation d, if d is not input.
#' @param baseline baseline rate parameter for random variables not given a value from Path
#' @param family string, either "bernoulli" or "poisson" is supported at this time, indicating the family of random variables to simulate in dataset
makeDataset <- function(Path,row.tree=NULL,col.tree=NULL,d=NULL,sd=100,baseline=0,family='bernoulli'){
  if ('treesim' %in% class(Path)){
    W <- Path$W
    V <- Path$V
    if (is.null(d)){
      d <- Path$D
    }
    Path <- Path$Paths
  } else {
    if (is.null(row.tree)){
      stop('if input Path is not a treesim object, must input row.tree')
    }
    if (is.null(col.tree)){
      stop('if input Path is not a treesim object, must input col.tree')
    }
    W <- treeBasis(row.tree)[,Path$row.node-ape::Ntip(row.tree),drop=F]
    V <- treeBasis(col.tree)[,Path$col.node-ape::Ntip(col.tree),drop=F]
  }
 
  if (is.null(d)){
    d <- rnorm(n=nrow(Path),sd=sd)
    d <- abs(d)*Path$orientation
    d <- diag(d)
  } else {
    if (is.null(dim(d))){
      if (length(d) != nrow(Path)){
        stop('if d is a vector, it must be of length equal to nrow(Path)')
      } else {
        if (!all(sign(d)==Path$orientation)){
          warning('sign(d) not all equal to Path$orientation - will be forced into equality')
          d <- abs(d)*Path$orientation
        }
        d <- diag(d)
      }
    }
  }
  m <- nrow(W)
  n <- nrow(V)
  X <- W %*% d %*% t(V)
  
  if (family=='bernoulli'){
    probs <- binomial(link='logit')$linkinv(X)
    probs[X==0] <- baseline
    X <- rbinom(m*n,1,c(probs)) %>% matrix(nrow=m,ncol=n,byrow=F)
  } else if (family=='poisson'){
    lambda <- poisson()$linkinv(X)
    lambda[X==0] <- baseline
    X <- rpois(m*n,c(lambda)) %>% matrix(nrow=m,ncol=n,byrow=F)    
  } else if (family=='gaussian'){
    
  } else {
    warning(paste('family',family,'not supported at this time. Outputting raw matrix, X'))
  }
  
  rownames(X) <- row.tree$tip.label
  colnames(X) <- col.tree$tip.label
  return(X)
}
