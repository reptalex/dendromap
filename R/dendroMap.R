
dendroMap <- function(X,row.tree,col.tree,forbidden.nodes=NULL,ncores=NULL,quantile.threshold=0.7){
  
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
  
  ### a function for internal use - find maximum nodeSeq
  seqstat <- function(Seq){
    if (is.null(Seq)){
      return(-Inf)
    } else {
      return(sum(abs(Seq$statistic)))
    }
  }
  
  row.nodemap <- dendromap:::makeNodeMap(row.tree)
  col.nodemap <- dendromap:::makeNodeMap(col.tree)
  W <- treeBasis(row.tree)
  V <- treeBasis(col.tree)
  U <- t(W) %*% X %*% V
  
  
  ###### Obtain null ecdf
  Unull <- t(W) %*% shuffleData(X) %*% V
  chisq <- log(unlist(Unull)^2)
  chisq <- chisq[chisq>-10]
  null.ecdf <- ecdf(chisq)
  rm(list=c('Unull','chisq','V','W'))
  gc()
  ######
  ###### we'll only consider node-pairs above a threshold - this threshold should be low
  threshold <- quantile(null.ecdf,quantile.threshold)
  threshold <- sqrt(exp(threshold)) ### only consider those with abs(U)>threshold
  U[abs(U)<threshold] <- 0
  ######
  
  done <- F
  iteration=0
  Seq <- NULL
  while (!done){
    iteration=iteration+1
    
    RowColSet <- dendromap:::makeRowColSet(U,forbidden.nodes)
    ######### Find node-seq maximizing objective
    if (is.null(ncores)){
      seqs <- lapply(RowColSet,dendromap:::findNodeSeq,U=U,
                     row.tree=row.tree,col.tree=col.tree,
                     row.nodemap=row.nodemap,col.nodemap=col.nodemap)
    } else {
      if (iteration==1){
        cl <- parallel::makeCluster(ncores)
        parallel::clusterEvalQ(cl,library(dendromap))
        parallel::clusterExport(cl,varlist=c('U','row.tree','col.tree',
                                             'row.nodemap','col.nodemap','getIndexSets'))

        seqs <- tryCatch(parallel::parLapply(cl,RowColSet,findNodeSeq),
                         error=function(e) e)
       if ('error' %in% class(seqs)){
          error.message <- paste('Error in findBestSeq. Iteration',iteration,
                                 'returned the following error:',as.character(seqs))
          parallel::stopCluster(cl)
          rm('cl')
          gc()
          stop(error.message)
        }
      }
    }
    max.nodeseq <- which.max(sapply(seqs,seqstat))
    if (length(max.nodeseq)==0 | is.na(max.nodeseq) | is.infinite(seqstat(seqs[[max.nodeseq]]))){
      stop(paste('max.nodeseq returned unacceptable value:',max.nodeseq))
    }
    
    ####### find ancestors
    Seq <- rbind(Seq,seqs[[max.nodeseq]])
    row.node <- Seq$row.node[1]
    col.node <- Seq$col.node[1]
    
    ancestors <- matchAncestors(U,row.node,col.node,row.tree,col.tree)
    if (!is.null(ancestors)){
      # ancestor.seq <- ...F(ancestors)
      # Seq <- rbind(Seq,)
    }
    ####### find descendants: 
    descendants <- matchDescendants(U,seqs[[max.nodeseq]],row.tree,col.tree,
                                    row.nodemap,col.nodemap,quantile.threshold)
    
  }
  
  which.max(sapply(seqs,seqstat))
}