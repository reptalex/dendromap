#' Match descendants using basalMaxima
matchDescendants <- function(U,Seq,row.tree,col.tree,
                             row.nodemap,col.nodemap,
                             quantile.threshold=0.95){
  orientation.filter <- function(contenders,Seq,row.nodemap,col.nodemap,start=1){
    for (i in start:nrow(Seq)){
      row.ix <- getIndexSets(Seq$row.node[i],row.nodemap)
      col.ix <- getIndexSets(Seq$col.node[i],col.nodemap)
      orientation <- sign(Seq$statistic[i])
      
      
      if (length(unlist(row.ix))==0 | length(unlist(col.ix))==0){
        terminal.row <- length(unlist(row.ix))==0 
        terminal.col <- length(unlist(col.ix))==0
        if (terminal.row & terminal.col){
          next
        } else if (terminal.row){
          ### remove this row.node from consideration of downstream column nodes
          contenders <- contenders[!(row.node==Seq$row.node[i] & col.node %in% unlist(col.ix))]
        } else {
          contenders <- contenders[!(col.node==Seq$col.node[i] & row.node %in% unlist(row.ix))]
        }
      } else {
        if (orientation==1){
          contenders <- contenders[((row.node %in% row.ix$pos) & (col.node %in% col.ix$pos))|
                                     ((row.node %in% row.ix$neg) & (col.node %in% col.ix$neg))]
        } else {
          contenders <- contenders[((row.node %in% row.ix$pos) & (col.node %in% col.ix$neg))|
                                     ((row.node %in% row.ix$neg) & (col.node %in% col.ix$pos))]
        }
      }
    }
    return(contenders)
  }
  
  ### this will involve basalMaxima for the descendants of our row.node
  ### Possibly with Metropolis-Hastings stochastic search in future versions
  m <- ape::Ntip(row.tree)
  n <- ape::Ntip(col.tree)
  row.ix <- getIndexSets(Seq$row.node[1],row.nodemap)
  col.ix <- getIndexSets(Seq$col.node[1],col.nodemap)
  
  f <- ecdf(c(abs(U)))
  u <- U[unlist(row.ix)-m,unlist(col.ix)-n,drop=F]
  u[f(abs(c(u)))<quantile.threshold] <- 0
  u[paste('node_',Seq$row.node[2],sep=''),
    paste('node_',Seq$col.node[2],sep='')] <- 0
  
  rows <- rownames(u)[which(u!=0) - (nrow(u))*(ceiling(which(u!=0)/nrow(u))-1)]
  cols <- colnames(u)[ceiling(which(u!=0)/nrow(u))]
  contenders <- data.table('row.node'=as.numeric(gsub('node_','',rows)),
                           'col.node'=as.numeric(gsub('node_','',cols)))
  ### set entries to 0 if they are on the (-) side of current nodes
  contenders <- orientation.filter(contenders,Seq,row.nodemap,col.nodemap)
  while (nrow(contenders)!=0){
    winner <- which.max(abs(U[cbind(contenders$row.node-m,contenders$col.node-n)]))
    Seq <- rbind(Seq,data.frame('row.node'=contenders$row.node[winner],
                                  'col.node'=contenders$col.node[winner],
                                  'statistic'=U[contenders$row.node[winner]-m,
                                                contenders$col.node[winner]-n]))
    contenders <- contenders[setdiff(1:nrow(contenders),winner)]
    if (nrow(contenders)!=0){
      contenders <- orientation.filter(contenders,Seq,row.nodemap,col.nodemap,
                                       start=nrow(Seq))
    }
  }
  return(Seq)
}
