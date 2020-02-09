devtools::install_github('reptalex/dendromap')
# devtools::install_github('reptalex/phylodt')
library(dendromap)


load('data/birds/bird_conts.Rdata')
load('data/birds/PTrees.Rdata')
row.tree <- PTrees[[1]]
col.tree <- readRDS('data/birds/continental_tree')

colnames(bird_conts)[c(2,5,6)] <- c('Eurasia','NAmerica','SAmerica')

# checking ----------------------------------------------------------------

X <- bird_conts[,2:6] %>% as.matrix
rownames(X) <- bird_conts$species
colnames(X) <- colnames(bird_conts)[2:6]

all(rownames(X) %in% row.tree$tip.label)
all(colnames(X) %in% col.tree$tip.label)


row.tree <- drop.tip(row.tree,setdiff(row.tree$tip.label,rownames(X)))

X <- X[row.tree$tip.label,col.tree$tip.label]


dm <- dendromap2(X,row.tree,col.tree,estimate_runtime = TRUE)
dm_comp <- dendromap2(X,row.tree,col.tree,estimate_runtime = TRUE,discard_contingency_zeros = TRUE)

saveRDS(dm,'data/birds/bird_dendromap')
saveRDS(dm_comp,'data/birds/bird_dendromap_nonzero_contingencies')
row.nodemap <- makeNodeMap(row.tree) #makes it easy to find pos/neg descendants for quick summaries
col.nodemap <- makeNodeMap(col.tree)
save(list=ls(),file='data/birds/bird_dendromap_workspace')


# summary -----------------------------------------------------------------


# Lin=3
# nds <- dm$Lineages[Lineage==Lin,row.node]
# spp <- vector(mode='list',length=length(nds))
# i=0
# for (nd in nds){
#   i=i+1
#   spp[[i]][[1]] <- phangorn::Descendants(row.tree,nd+1,'tips')[[1]] %>% row.tree$tip.label[.] %>% sort
#   spp[[i]][[2]] <-  phangorn::Descendants(row.tree,nd+row.nodemap[node==nd,pos]+1,'tips')[[1]] %>% row.tree$tip.label[.] %>% sort
# }
# 
# 
# i=1
# continents <- phangorn::Descendants(col.tree,dm$Lineages[Lineage==Lin][i,col.node],'tips')[[1]] %>% col.tree$tip.label[.]
# X[unlist(spp[[i]]),continents]




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
    
    # if (length(grp1)==0){
    #   grp1 <- phangorn::Descendants(col.tree,colnds[k],'tips')[[1]][1] %>% dm$row.tree$tip.label[.]
    # }
    # if (length(grp2)==0){
    #   grp2 <- phangorn::Descendants(col.tree,colnds[k],'tips')[[1]][1] %>% dm$row.tree$tip.label[.]
    # }
    
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
    
    
    base::cat(paste('\n\n\n ------------------------------------------------------------ \n event',k,'for lineage',i,' \n Row tree species split - spp1|spp2'))
    base::cat(paste('\n',rowsplit))
    base::cat(paste('\n',length(spp1),'species in spp1 |',length(spp2),'species in spp2'))
    base::cat(paste('\n\n Column tree split grp1|grp2'))
    base::cat(paste('\n',colsplit))
    base::cat('\n\n\n ------ Contingency Table ------\n\n     grp1     grp2\nspp1  ',signif(tbl[1,1],2),' | ',signif(tbl[1,2],2))
    base::cat('\n-------------------------')
    base::cat('\nspp2  ',signif(tbl[2,1],2),' | ',signif(tbl[2,2],2))
    
    tables[[k]] <- tbl
    tiplabels[[k]] <- list('species1'=spp1,'species2'=spp2,
                           'group1'=grp1,'group2'=grp2)
  }
  return(list('tables'=tables,
              'tiplabels'=tiplabels))
}



summarize_lineage(1,dm,row.nodemap,col.nodemap)
summarize_lineage(2,dm,row.nodemap,col.nodemap)


summarize_lineage(5,dm,row.nodemap,col.nodemap)

for (i in 1:max(dm$Lineages$Lineage)){
  summarize_lineage(i,dm,row.nodemap,col.nodemap)
}



for (i in 1:max(dm_comp$Lineages$Lineage)){
  summarize_lineage(i,dm_comp,row.nodemap,col.nodemap)
}
