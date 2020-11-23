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



# filter to hominoidea ----------------------------------------------------



dm <- dendromap(Data,row.tree,col.tree,ncores=7)  

