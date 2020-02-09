#' Summarize lineage distribution in dendromap object
#' @export
#' @param dm \code{\link{dendromap}} object
#' @param row.nodemap optional \code{row.tree} nodemap from \code{\link{makeNodeMap}}
#' @param col.nodemap optional \code{col.tree} nodemap from \code{\link{makeNodeMap}}
summarize_lineage <- function(i,dm,row.nodemap=NULL,col.nodemap=NULL){
  if (is.null(row.nodemap)){
    row.nodemap <- makeNodeMap(dm$row.tree)
  }
  if (is.null(col.nodemap)){
    col.nodemap <- makeNodeMap(dm$col.tree)
  }
  Lineage <- dm$Lineages[Lineage==i]
  rownds <- Lineage[,row.node]
  colnds <- Lineage[,col.node]
  
  nr <- ape::Ntip(row.tree)
  nc <- ape::Ntip(col.tree)
  
  tables <- vector(mode='list',length=nrow(Lineage))
  tiplabels <- vector(mode='list',length=nrow(Lineage))
  for (k in 1:nrow(Lineage)){
    
    rchildren <- row.tree$edge[row.tree$edge[,1]==rownds[k],2]
    if (rchildren[1]>nr){
      spp1 <- phangorn::Descendants(row.tree,rownds[k]+1,'tips')[[1]] %>% dm$row.tree$tip.label[.]
    } else {
      spp1 <- dm$row.tree$tip.label[rchildren[1]]
    }
    if (rchildren[2]>nr){
      spp2 <-  phangorn::Descendants(row.tree,rownds[k]+row.nodemap[node==rownds[k],pos]+1,'tips')[[1]] %>% dm$row.tree$tip.label[.]
    } else {
      spp2 <- dm$row.tree$tip.label[rchildren[2]]
    }
    
    children <- col.tree$edge[col.tree$edge[,1]==colnds[k],2]
    if (children[1]>nc){
      grp1 <- phangorn::Descendants(col.tree,colnds[k]+1,'tips')[[1]] %>% dm$col.tree$tip.label[.]
    } else {
      grp1 <- dm$col.tree$tip.label[children[1]]
    }
    if (children[2]>nc){
      grp2 <- phangorn::Descendants(col.tree,colnds[k]+col.nodemap[node==colnds[k],pos]+1,'tips')[[1]] %>% dm$col.tree$tip.label[.]
    } else {
      grp2 <- dm$col.tree$tip.label[children[2]]
    }
    
    tbl <- matrix(NA,nrow=2,ncol=2)
    rownames(tbl) <- c('spp1','spp2')
    colnames(tbl) <- c('grp1','grp2')
    tbl[1,1] <- mean(dm$Data[spp1,grp1],na.rm=T)
    tbl[1,2] <- mean(dm$Data[spp1,grp2],na.rm=T)
    tbl[2,1] <- mean(dm$Data[spp2,grp1],na.rm=T)
    tbl[2,2] <- mean(dm$Data[spp2,grp2],na.rm=T)
    rowsplit=paste(paste(unique(sample(spp1,2,replace=T)),collapse=','),'|',
                   paste(unique(sample(spp2,2,replace=T)),collapse=','))
    colsplit=paste(paste(grp1,collapse=','),'|', paste(grp2,collapse=','))
    
    
    base::cat(paste('\n\n\n ------------------------------------------------------------ \n event',k,'for lineage',i))
    base::cat('\nP=',signif(Lineage$P[k],4))
    base::cat(paste('\n----------------------------------\n\nRow tree species split - species1|species2 \n',rowsplit))
    base::cat(paste('\n',length(spp1),'species in species1 |',length(spp2),'species in species2'))
    base::cat(paste('\n----------------------------------\n\nColumn tree split - group1|group2 \n',colsplit))
    base::cat('\n\n\n ------ Contingency Table ------\n\n     grp1     grp2\nspp1  ',signif(tbl[1,1],2),' | ',signif(tbl[1,2],2))
    base::cat('\n-------------------')
    base::cat('\nspp2  ',signif(tbl[2,1],2),' | ',signif(tbl[2,2],2))
    
    tables[[k]] <- tbl
    tiplabels[[k]] <- list('species1'=spp1,'species2'=spp2,
                           'group1'=grp1,'group2'=grp2)
  }
  return(list('tables'=tables,
              'tiplabels'=tiplabels))
}