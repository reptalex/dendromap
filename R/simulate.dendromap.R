simulate.dendromap <- function(dm,...){
  m <- length(dm$row.tree$tip.label)
  n <- length(dm$col.tree$tip.label)
  
  Lineages <- dm$Lineages
  setkey(Lineages,lineage_id,row.edge)
  
  # X <- Matrix::spMatrix(m,n)
  X <- matrix(0,m,n)
  rownames(X) <- dm$row.tree$tip.label
  colnames(X) <- dm$col.tree$tip.label
  ids <- Lineages$lineage_id %>% unique
  for (id in ids){
    Ln <- Lineages[lineage_id==id]
    n_grains <- nrow(Ln)
    for (ll in 1:n_grains){
      
      re <- Ln$row.edge[ll]
      ce <- Ln$col.edge[ll]
      is <- dm$rowEdgeTips[edge==re,list(seq(min,max))]$V1
      js <- dm$colEdgeTips[edge==ce,list(seq(min,max))]$V1
      if (ll==1){
        X[is,js] <- 1
        rowset <- is
        colset <- js
      } else {
        comp_set <- setdiff(colset,js)
        X[is,comp_set] <- 0
      }
    }
  }
  return(X)
}
