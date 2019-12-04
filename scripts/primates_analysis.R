library(dendromap)
library(ggplot2)
library(biomformat)
library(Biostrings)
library(ggpubr)
library(tictoc)

##### BasalMaxima analysis of primate data

B <- read_biom('data/primates/final.withtax.min10.biom')
Data <- biom_data(B)
taxonomy <- data.frame('sequence'=sapply(B$rows,getElement,'id'),
                       'taxonomy'=sapply(B$rows,FUN=function(x) paste(x$metadata$taxonomy,collapse='; ')))
taxonomy$taxonomy <- as.character(taxonomy$taxonomy)
taxonomy <- as.data.table(taxonomy)
setkey(taxonomy,sequence)
X <- read.csv('data/primates/Reduced_mapping_4persp_nohyb_RC2036.txt',sep = '\t',stringsAsFactors = F)
names(X)[1] <- 'SampleID'
X <- as.data.table(X)
setkey(X,SampleID)
Data <- Data[,match(X$SampleID,colnames(Data))]
all.equal(colnames(Data),X$SampleID)



# trees -------------------------------------------------------------------
seqs <- as.list(rownames(Data))
paste(paste('>',seqs,sep=''),seqs,sep='\n') %>%
  write(file='data/primates/final_withtax_min10.fasta',sep='\n')
seqs <- Biostrings::DNAStringSetList(seqs)
write.dna(seqs,file='data/primates/final_withtax_min10.fasta')

col.tree <- read.tree('data/primates/furcated_host_tree.tre')
row.tree <- read.tree('data/primates/final_withtax_min10.tre')
all(row.tree$tip.label %in% rownames(Data))
all(col.tree$tip.label %in% colnames(Data))
Data <- Data[row.tree$tip.label,]
Data <- Data[,col.tree$tip.label]



# Is there signal? --------------------------------------------------------
W <- treeBasis(row.tree)
V <- treeBasis(col.tree)

U <- t(W) %*% Data %*% V

shuffleData <- function(Data,rows=TRUE,cols=TRUE){
  if (rows){
    ix.rows <- sample(nrow(Data))
  } else {
    ix.rows <- 1:nrow(Data)
  }
  if (cols){
    ix.cols <- sample(ncol(Data))
  } else {
    ix.cols <- 1:ncol(Data)
  }
  return(Data[ix.rows,ix.cols])
}

Unull <- t(W) %*% shuffleData(Data) %*% V


DF <- data.table('z'=c(log(c(as.matrix(U))^2),
                       log(c(as.matrix(Unull))^2)),
                 'Dataset'=rep(c('Primates','Null'),each=nrow(U)*ncol(U)))

ggplot(DF[z>-10],aes(z,by=Dataset,fill=Dataset))+
  geom_histogram(bins=100,alpha=0.4,position='identity')

z <- DF[Dataset=='Primates' & z>-10,z]
znull <- DF[Dataset=='Null' & z>-10,z]



# dendromap ---------------------------------------------------------------


tic()
dm <- dendromap(Data,row.tree,col.tree,W=W,V=V,Pval_threshold = 0.001,ncores=7)  
toc()
#24.75s
# saveRDS(dm,'data/primates/dendromap_P8e-4_threshold')
saveRDS(dm,'data/primates/dendromap_P1e-3_threshold')
# 
# tic()
# dm <- dendromap(Data,row.tree,col.tree,W=W,V=V,Pval_threshold = 0.005,ncores=7)  
# toc() ## here, the max clique step took most of the time.
# 1086 RCs had P<=Pval_threshold at Pval_threshold=0.005
# RCmap has 4160 rows
# There are 666 terminal nodes. Traversing RC tree from all terminal nodes.
# Found 2264 RC sequences for 2561716 pairs. Checking joinability of pairs.
# 146.14s
# saveRDS(dm,'data/primates/dendromap_P5e-3_threshold')

Lineages <- dm$Lineages
Lineages[,n:=.N,by=Lineage]

Lineages[n==max(n),min(row.node)] ### this is the descendant node of our biggest lineage
setkey(Lineages,n,Lineage)

ll=2
row.descendants <- row.tree$tip.label[phangorn::Descendants(row.tree,Lineages[Lineage==ll,min(row.node)],'tips')[[1]]]
col.descendants <- col.tree$tip.label[phangorn::Descendants(col.tree,Lineages[Lineage==ll,min(col.node)],'tips')[[1]]]
taxonomy[sequence %in% row.descendants]
# X[SampleID %in% col.descendants,c('Species_tree','Diet','Phyl_Group')]


taxonomic.summary <- function(dm,taxonomy,lineage=1){
  row.descendants <- row.tree$tip.label[phangorn::Descendants(row.tree,dm$Lineages[Lineage==lineage,min(row.node)],'tips')[[1]]]
  ix <- as.character(unlist(taxonomy[,1])) %in% row.descendants
  tx <- taxonomy[which(ix),2] %>% unlist
  tx2 <- taxonomy[which(!ix),2] %>% unlist
  pp <- phylofactor::uniqueTaxa(tx,tx2)
  return(unique(pp))
}

plot(dm)

plot_lineage <- function(dm,lineage=1,y=NULL,trans='none',row.tip.label=FALSE,col.tip.label=TRUE,...){
  if (is.null(y)){
    if (is.null(dm$Data)){
      ## plot just the lineages on row and column trees, with tip labels
    } else {
      ## use dm$Data for a plot of subset containing lineage
      spp <- phangorn::Descendants(dm$row.tree,dm$Lineages[Lineage==lineage,min(row.node)],'tips')[[1]]
      spp <- dm$row.tree$tip.label[spp]
      
      xx <- dm$Data[spp,]
      if (trans!='none'){
        if (trans=='log'){
          min.val <- tryCatch(min(xx[xx!=0]),error=function(e) e)
          if (min.val<0 | 'error' %in% class(min.val)){
            stop(paste('error produced a finding min(dm$Data[dm$Data!=0]): min.val=',e))
          } else {
            xx[xx==0] <- 0.65*min.val
            xx <- log(xx)
          }
        } else {
          stop('trans must be either "none" or "log"')
        }
      }
    }
  } else {
    spp <- phangorn::Descendants(dm$row.tree,dm$Lineages[Lineage==lineage,min(row.node)],'tips')[[1]]
    spp <- dm$row.tree$tip.label[spp]
    xx <- y[spp,]
  }
  if (exists('xx')){
    dd <- dm[setdiff(names(dm),'Data')]
    class(dd) <- 'dendromap'
    dd$Data <- xx
    dd$Lineages <- dd$Lineages[Lineage==lineage]
    groups <- lapply(dd$Lineages$row.node,FUN=function(x,tr) tr$tip.label[phangorn::Descendants(tr,x,'tips')[[1]]],tr=dd$row.tree)
    dd$row.tree <- ape::drop.tip(dd$row.tree,setdiff(dd$row.tree$tip.label,spp))
    ## since we trimmed the row tree, we need to go back and replace the row.nodes
    ## to match the mrca's
    
    dd$Lineages$row.node <- sapply(groups,FUN=function(x,tr) ape::getMRCA(tr,x),tr=dd$row.tree)
    
    
    col.root <- dd$Lineages[,min(col.node)]
    col.tips <- phangorn::Descendants(dd$col.tree,col.root,'tips')[[1]]
    
    dd$Data <- dd$Data[,col.tips]
    plot(dd,...)
  }
}


plot_lineage(dm,1,trans='log',col.tr.left=0.48,col.tr.width=0.49)

plot_lineage(dm,2,trans='log',col.tr.left=0.48,col.tr.width=0.49)

taxonomic.summary(dm,taxonomy)
