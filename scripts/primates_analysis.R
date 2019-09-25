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
dm <- dendromap(Data,row.tree,col.tree,W=W,V=V,Pval_threshold = 0.0008)  
toc()

Lineages <- dm$Lineages
Lineages[,n:=.N,by=Lineage]

Lineages[n==max(n),min(row.node)] ### this is the descendant node of our biggest lineage
setkey(Lineages,n,Lineage)

ll=4
row.descendants <- row.tree$tip.label[phangorn::Descendants(row.tree,Lineages[Lineage==ll,min(row.node)],'tips')[[1]]]
col.descendants <- col.tree$tip.label[phangorn::Descendants(col.tree,Lineages[Lineage==ll,min(col.node)],'tips')[[1]]]
taxonomy[sequence %in% row.descendants]
X[SampleID %in% col.descendants,c('Species_tree','Diet','Phyl_Group')]
