#' Identify rc pairs for which we have data to compare significance
#' @export
#' @param N dataset
#' @param W row.tree treeBasis
#' @param V col.tree treeBasis
#' @param cl cluster
#' @param verbose logical, whether or not to display message (used inside dendromap for more useful updates of algorithm progress)
#' @examples
#' 
#' # If a row tree has descendant clades represented by {clade1,clade2} and col tree has descendant clades {continent1,continent2}
#' # Then a row-col node pair is considered 'incomparable' if any row or column of its contingency table has all zeros.
#' # e.g.: incomparable:
#' #       cnt1   cnt2  |        cnt1   cnt2
#' # spp1   0      0    |  spp1   0     0.5
#' # spp2   1     0.5   |  spp2   0      1
#' #
#' #e.g. comparable:
#' #
#' #       cnt1   cnt2  |        cnt1   cnt2
#' # spp1   0      1    |  spp1   0     0.5
#' # spp2   1     0.5   |  spp2   1      1
idComparables <- function(N,W,V,cl=NULL,verbose=TRUE){
  if (verbose){
    base::cat('Identifying comparable species|species x continent|continent splits')
  }
  isComparable <- function(i,j,N.=N,W.=W,V.=V){
    w <- W[,i]
    v <- V[,j]
    spp1 <- which(w>0)
    spp2 <- which(w<0)
    grp1 <- which(v>0)
    grp2 <- which(v<0)
    spp_present <- any(N[spp1,c(grp1,grp2)]!=0) & any(N[spp2,c(grp1,grp2)]!=0) #false if either species is missing
    sites_present <- any(N[c(spp1,spp2),grp1]!=0) & any(N[c(spp1,spp2),grp2]!=0)
    if (spp_present & sites_present){
      return(TRUE)
    } else {
      return(FALSE)
    }
  }
  
  n <- ncol(W)
  m <- ncol(V)
  tbl <- data.table('i'=rep(1:n,times=m),
                    'j'=rep(1:m,each=n))
  if (is.null(cl)){
    comparable=apply(tbl,1,FUN=function(x) isComparable(x[1],x[2]))
  } else {
    parallel::clusterExport(cl,varlist=c('N','W','V'),envir = environment())
    parallel::clusterExport(cl,varlist='isComparable')
    comparable=parallel::parApply(cl,tbl,1,FUN=function(x) isComparable(x[1],x[2]))
  }
  tbl[,comparable:=comparable]
  return(tbl)
}
