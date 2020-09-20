#' internal function to find joinable RC sequences
#' @export
#' @param Seqs output from \code{\link{rc_seqs}}
#' @param rc_table output from \code{\link{makeRCtable}}
#' @param Row_Descendants named \code{getIndexSets} of all row.nodes
#' @param Col_Descendants named \code{getIndexSets} of all col.nodes
#' @param method either 'node' or 'edge
#' @param cl optional cluster with dendromap loaded
find_joinables <- function(Seqs,rc_table,Row_Descendants,Col_Descendants,method='node',cl=NULL){
  setkey(rc_table,rc_index)
  n <- length(Seqs)
  tbl <- data.table('seq1'=rep(1:(n-1),times=(n-1):1),key='seq1')
  tbl[,seq2:=(seq1+1):n,by=seq1]
  # base::cat(paste('\nFound',length(Seqs),'RC sequences for',nrow(tbl),'pairs. Checking joinability of pairs.'))
  getBranchPoint <- function(seq1,seq2,Seqs){
    intrsct <- intersect(Seqs[[seq1]],Seqs[[seq2]])
    branch_points1 <- !Seqs[[seq1]]%in%intrsct
    branch_points2 <- !Seqs[[seq2]]%in%intrsct
    
    ix1 <- Seqs[[seq1]][min(which(branch_points1))]
    ix2 <- Seqs[[seq2]][min(which(branch_points2))]
    return(c('ix1'=ix1,'ix2'=ix2))
  }
  rc_ix <- unique(unlist(Seqs))
  rctbl <- rc_table[rc_index %in% rc_ix]
  
  ### trim table
  A <- sapply(Seqs,FUN=function(a,b) as.numeric(b %in% a),b=rc_ix) %>%
    Matrix::Matrix(sparse=T)# will have one colum for every Seq
  ix <- (Matrix::t(A) %*% A)[as.matrix(tbl)]>0 ### here we can gpu-compute t(A) %*% A - VERY useful for large tbl & sparse mat
  rm('A')
  tbl <- tbl[ix]
  rm('ix')
  if (nrow(tbl)>0){
    branch_points <- mapply(FUN=getBranchPoint,seq1=tbl$seq1,seq2=tbl$seq2,
                            MoreArgs = list('Seqs'=Seqs)) %>% t %>% as.data.table()
    branch_points[,ix:=1:.N]
    colnames(branch_points)[1] <- 'rc_index'
    setkey(branch_points,rc_index)
    branch_points <- rctbl[,c('rc_index','row.node','col.node')][branch_points]
    names(branch_points)[1:4] <- c('ix1','rn1','cn1','rc_index')
    setkey(branch_points,rc_index)
    branch_points <- rctbl[,c('rc_index','row.node','col.node')][branch_points]
    names(branch_points)[1:3] <- c('ix2','rn2','cn2')
    if (is.null(cl)){
      joinability <- apply(as.matrix(branch_points[,c('rn1','rn2','cn1','cn2')]),1,
                           FUN=function(x,r,c) check_joinable(x[1],x[2],x[3],x[4],r,c),
                           r=Row_Descendants,c=Col_Descendants)
    } else {
      joinability <- parallel::parApply(cl=cl,X=as.matrix(branch_points[,c('rn1','rn2','cn1','cn2')]),1,
                                         FUN=function(x,r,c) check_joinable(x[1],x[2],x[3],x[4],r,c),
                                         r=Row_Descendants,c=Col_Descendants)
    }
    joinability <- joinability[order(branch_points$ix)]
    tbl[,joinability:=joinability]
    setkey(rc_table,row.node,col.node)
  } else {
    tbl <- NULL
  }
  return(tbl)
}
