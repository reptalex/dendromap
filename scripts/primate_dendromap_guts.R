library(dendromap)
library(ggplot2)
library(biomformat)
library(Biostrings)
library(ggpubr)
library(profvis)
library(tictoc)

##### BasalMaxima analysis of primate data

B <- read_biom('data/primates/final.withtax.min10.biom')
X <- biom_data(B)
taxonomy <- data.frame('sequence'=sapply(B$rows,getElement,'id'),
                       'taxonomy'=sapply(B$rows,FUN=function(x) paste(x$metadata$taxonomy,collapse='; ')))
taxonomy$taxonomy <- as.character(taxonomy$taxonomy)
taxonomy <- as.data.table(taxonomy)
setkey(taxonomy,sequence)
M <- read.csv('data/primates/Reduced_mapping_4persp_nohyb_RC2036.txt',sep = '\t',stringsAsFactors = F)
names(M)[1] <- 'SampleID'
M <- as.data.table(M)
setkey(M,SampleID)
X <- X[,match(M$SampleID,colnames(X))]
all.equal(colnames(X),M$SampleID)



# trees -------------------------------------------------------------------
# seqs <- as.list(rownames(X))
# paste(paste('>',seqs,sep=''),seqs,sep='\n') %>%
  # write(file='data/primates/final_withtax_min10.fasta',sep='\n')
# seqs <- Biostrings::DNAStringSetList(seqs)
# write.dna(seqs,file='data/primates/final_withtax_min10.fasta')

col.tree <- read.tree('data/primates/furcated_host_tree.tre')
row.tree <- read.tree('data/primates/final_withtax_min10.tre')
all(row.tree$tip.label %in% rownames(X))
all(col.tree$tip.label %in% colnames(X))
X <- X[row.tree$tip.label,]
X <- X[,col.tree$tip.label]



# Is there signal? --------------------------------------------------------
W <- treeBasis(row.tree)
V <- treeBasis(col.tree)


# dendromap ---------------------------------------------------------------
set.seed(1)
ncores=7
Pval_threshold=0.005
n_sim=NULL
cl <- parallel::makeCluster(ncores)
parallel::clusterEvalQ(cl,library(dendromap))

row.nodemap <- dendromap:::makeNodeMap(row.tree)
col.nodemap <- dendromap:::makeNodeMap(col.tree)

rc_table <- makeRCtable(X,row.tree,col.tree,W,V,n_sim)
rc_table <- rc_table[P<max(P)]

rc_table <- rc_table[P<=Pval_threshold]
# 1000 RCs had P<=Pval_threshold at Pval_threshold=0.005

row.nodes <- unique(rc_table$row.node)
col.nodes <- unique(rc_table$col.node)
Row_Descendants <- lapply(row.nodes,getIndexSets,row.nodemap) %>%
  lapply(FUN=function(x,a) lapply(x,intersect,a),a=row.nodes)
Col_Descendants <- lapply(col.nodes,getIndexSets,col.nodemap) %>%
  lapply(FUN=function(x,a) lapply(x,intersect,a),a=col.nodes)

names(Row_Descendants) <- row.nodes
names(Col_Descendants) <- col.nodes
RCmap <- makeRCMap(rc_table,Row_Descendants,Col_Descendants)




# find lineages - where things get NP hard --------------------------------
nds <- unique(RCmap[terminal==TRUE,descendant])
base::cat('\nThere are',length(nds),'terminal nodes. Traversing RC tree from all terminal nodes.')

#### RC SEQUENCES = ANT PATHS
Seqs <- parallel::parLapply(cl=cl,nds,rc_seqs,RCmap) %>% unlist(recursive=FALSE) %>% unique

#### JOINABILITY
tbl <- dendromap:::find_joinables(Seqs,rc_table,Row_Descendants,Col_Descendants,cl)
# Found 2006 RC sequences for 2011015 pairs. Checking joinability of pairs.
nrow(tbl) # 175171 - 2004 RC seqs with 175,171 edges

#### JOINABILITY GRAPH
seq.indexes <- unique(unlist(tbl[,c('seq1','seq2')]))
nverts <- length(seq.indexes) # 2004
VertMap <- data.table('seq'=seq.indexes,'vert'=paste('v',seq.indexes,sep='_'),key='vert')

joinables <- tbl[joinability==TRUE]
nedge <- nrow(joinables)
joinable.edges <- split(joinables[,c('seq1','seq2')],seq(nrow(joinables))) %>% unlist
joinable.edges <- VertMap[match(joinable.edges,seq),vert]

base::cat(paste('\nFound',nedge,'joinable RC sequences containing',nverts,' joinable pairs. \n...how hard is it to find maximal cliques in a graph of',length(unique(joinable.edges)),'vertices and',nedge,'edges?'))
G <- igraph::make_graph(joinable.edges,directed = F)

# save(list=ls(),file='data/primates/primate_dendromap_guts_workspace')


# Graph analysis ----------------------------------------------------------

#### This graph is too large to solve in a human time - speeding up this max clique problem will be necessary for our package
load('data/primates/primate_dendromap_guts_workspace')

plot(G,pch=16)
### G has many subgraphs

sG <- igraph::clusters(G)

sG.vertices <- lapply(1:sG$no,FUN=function(a,m) names(which(m==a)),m=sG$membership)

SubGraphs <- lapply(sG.vertices,FUN=function(v,graph) igraph::induced_subgraph(G,v),graph=G)
subgraph_sizes <- sapply(sG.vertices,length)
lns <- sapply(sG.vertices,length)
g <- SubGraphs[[which.max(lns)]]
g2 <- SubGraphs[[order(lns,decreasing = T)[2]]]
g3 <- SubGraphs[[order(lns,decreasing = T)[3]]]
### the biggest subgraph has 323 vertices and 21,425 edges
par(mfrow=c(1,3))
plot(g,vertex.size=5,vertex.label=NA,
     vertex.shape='sphere',vertex.color='black',edge.color='steelblue')
plot(g2,vertex.size=5,vertex.label=NA,
     vertex.shape='sphere',vertex.color='black',edge.color='steelblue')
plot(g3,vertex.size=5,vertex.label=NA,
     vertex.shape='sphere',vertex.color='black',edge.color='steelblue')

saveRDS(g,'data/primates/big_graph')
saveRDS(g2,'data/primates/big_graph2')
saveRDS(g3,'data/primates/big_graph3')




# iterative joining of compatible sisters ---------------------------------

### two sequences 
### ...-C-D-A1-A2-...-An
### ...-C-D-B1-B2-...-Bm
### are compatible sequences iff (1) An and Bm share no common descendants.
### and (2) the sequences are equal for the first k entries and unequal afterwards.
roots <- function(g,Seqs.=Seqs,VertMap.=VertMap){
  VM <- VertMap[match(igraph::vertex_attr(g)$name,vert)]
  names(Seqs) <- 1:length(Seqs)
  seqs <- Seqs[VM$seq]
  sapply(seqs,getElement,1) %>% unique %>% return
}

lapply(SubGraphs,roots) %>% sapply(length) ## all 1: every subgraph is defined by a common root
sapply(SubGraphs,roots) %>% unique %>% length ## all subgraph roots are distinct.
## are any subgraphs' roots anc-dec of another on microbe tree?



# RC seq tree -------------------------------------------------------------

### we need to make a tree for our seqs. They radiate from the root
### Each possibility - unique disagremeent in sequence - will be a node.
VM <- VertMap[match(igraph::vertex_attr(g)$name,vert)]
names(Seqs) <- paste('v_',1:length(Seqs),sep='')
seqs <- Seqs[VM$seq] #we'll whittle this down as we go.
n <- length(seqs) ## number of tips in our seq-tree

###### MAKE RC-SEQ TREE
lns <- sapply(seqs,length)
k <- max(lns)
pad <- function(s,k.=k) paste(c(s,rep('none',(k-length(s)))),collapse=';')
tr <- phylofactor::taxaTree(sapply(seqs,pad))
tr$tip.label <- names(tr$tip.label)
## the names of tip labels now contain the sequence name
## and the tip labels themselves contain the rc_index sequences
nnds <- ape::Nnode(tr)

setkey(rc_table,rc_index)
rcs <- rc_table[rc_index %in% unique(unlist(seqs))]
ix <- as.character(rcs$rc_index)
seq_scores <- -log(rcs$P)
names(seq_scores) <- ix
compute.score <- function(clq,seqs,seq_scores){
  nms <- unique(unlist(seqs[names(clq)]))
  sum(seq_scores[as.character(nms)]) %>% return
}



cl <- makeCluster(7)
clusterEvalQ(cl,library(dendromap))
tic()
losers <- NULL
for (nd in (n+nnds):(n+1)){
  kids <- tr2$edge[tr2$edge[,1]==nd,2]
  ss <- tr2$tip.label[phangorn::Descendants(tr2,nd,type = 'tips')[[1]]]
  n1=length(ss)
  ss <- setdiff(ss,losers)
  base::cat(paste('\nnode=',nd,' n1=',n1,' n=',length(ss),sep=''))
  
  subg <- igraph::induced_subgraph(G,ss)
  cliques <- igraph::max_cliques(subg)
  base::cat('\n---finished_cliques')
  scores <- parSapply(cl,cliques,compute.score,seqs,seq_scores)
  winner <- scores==max(scores)
  losers <- c(losers,unique(unlist(sapply(cliques[!winner],names))))
}
toc()

### winner: 
### cliques[winner]
# + 8/95 vertices, named, from 0efa86a:
# [1] v_1272 v_471  v_1549 v_314  v_958  v_1513 v_1981 v_1542
### score:
# [1] 92.54756

rclique <- function(g,verts=NULL,scores=NULL){
  if (is.null(verts)){
    verts <- igraph::vertex_attr(g)$name
  }
  if (is.null(scores)){
    scores <- rep(1,length(verts))
    names(scores) <- verts
  }
  A <- igraph::as_adjacency_matrix(g)
  nds <- sample(verts,1,prob=scores)
  remaining <- names(which(A[nds,]>0))
  while (length(remaining)>0){
    nds <- c(nds,sample(remaining,1,prob=scores[match(remaining,verts)]))
    connections <- (matrix(1,ncol=length(nds)) %*% A[nds,])[1,]
    remaining <- names(which(connections==length(nds)))
  }
  return(nds)
}

MH_clique <- function(g,nreps,seq_scores,lambda=1){
  verts <- igraph::vertex_attr(g)$name
  nd_list <- lapply(1:nreps,FUN=function(n,g,v,s) rclique(g,v,s),
                    g=g,v=verts,s=seq_scores[verts]^lambda)
  scores <- sapply(nd_list,FUN=function(clq,p) sum(seq_scores[clq]))
  winner <- which.max(scores)
  return(list('clique'=nd_list[[winner]],'score'=scores[[winner]]))
}


## REQUIRES:
# seq_scores with vertex names
# VertMap
# rclique exported
cl <- parallel::makeCluster(7)
parallel::clusterEvalQ(cl,library(dendromap))
parallel::clusterExport(cl,'rclique')
graph_size_threshold=30
big_graphs <- SubGraphs[subgraph_sizes>graph_size_threshold]
small_graphs <- SubGraphs[subgraph_sizes<=graph_size_threshold]
seq_scores <- parSapply(cl,Seqs,FUN=function(s,rc_table) -sum(log(rc_table[rc_index %in% s,P])),rc_table)
names(seq_scores) <- names(Seqs)
tic()
Cliques <- parSapply(cl,big_graphs,MH_clique,
                     nreps=100,seq_scores)
small_cliques <- sapply(small_graphs,igraph::max_cliques)
toc()
## 2.09s
parallel::stopCluster(cl)
rm('cl')


names(Seqs) <- paste('v_',1:length(Seqs),sep='')
rc_scores <- rc_table[,-log(P)]
names(rc_scores) <- rc_table$rc_index

max_clique_SA <- function(g,Seqs,rc_scores,seq_scores,time.limit=1,nreps=10){
  
  rclique <- function(g,verts,scores){
    if (is.null(verts)){
      verts <- igraph::vertex_attr(g)$name
    }
    if (is.null(scores)){
      scores <- rep(1,length(verts))
      names(scores) <- verts
    }
    A <- igraph::as_adjacency_matrix(g)
    nds <- sample(verts,1,prob=scores)
    remaining <- names(which(A[nds,]>0))
    while (length(remaining)>0){
      nds <- c(nds,sample(remaining,1,prob=scores[match(remaining,verts)]))
      connections <- (matrix(1,ncol=length(nds)) %*% A[nds,])[1,]
      remaining <- names(which(connections==length(nds)))
    }
    return(nds)
  }
  getScore <- function(s,Seqs,rc_scores)    sum(rc_scores[as.character(unique(unlist(Seqs[s])))])
  MH_clique <- function(g,Seqs,nreps,rc_scores,seq_scores,lambda=1){
    verts <- igraph::vertex_attr(g)$name
    nd_list <- lapply(1:nreps,FUN=function(n,g,v,s) rclique(g,v,s),
                      g=g,v=verts,s=seq_scores[verts]^lambda)
    scores <- sapply(nd_list,getScore,Seqs,rc_scores)
    winner <- which.max(scores)
    return(list('clique'=nd_list[[winner]],'score'=scores[[winner]]))
  }
  verts <- igraph::vertex_attr(g)$name
  A <- igraph::as_adjacency_matrix(g)
  one <- matrix(1,ncol=ncol(A))
  
  ### initialize with MH_clique
  clq <- MH_clique(g,Seqs,nreps,rc_scores,seq_scores)
  nds <- clq$clique
  start.time <- Sys.time()
  runtime <- Sys.time()-start.time
  
  i=0
  maxScore <- clq$score
  maxClique <- clq$clique
  while (runtime<=time.limit){
    remaining <- setdiff(verts,nds)
    neighbors <- matrix(1,ncol=length(nds)) %*% A[nds,]
    candidates <- neighbors[,remaining]
    candidates <- candidates[candidates>=(length(nds)-1)]
    if (length(candidates)==0){
      break
    } else if (length(candidates)>1){
      selected.candidate <- sample(candidates,1,prob=seq_scores[names(candidates)])
    } else {
      selected.candidate <- candidates
    }
    candidate.degree <- as.numeric(selected.candidate)
    selected.candidate <- names(selected.candidate)
    selected.candidate.score <- seq_scores[selected.candidate]
    if (candidate.degree<length(nds)){ 
      #some node(s) are not neighbors with our candidate - these are candidates for removal
      removal <- names(which(A[nds,selected.candidate]==0))
      if (length(removal)>1){
        removal <- sample(removal,1,prob=seq_scores[removal])
      }
    } else {
      removal <- sample(nds,1,prob=seq_scores[nds])
    }
    
    temperature <- max(as.numeric(time.limit-runtime),0)
    acceptance.probability <- exp(-(seq_scores[removal]-seq_scores[selected.candidate])/temperature)
    swap <- as.logical(runif(1)<acceptance.probability)
    if (swap){
      nds[nds==removal] <- selected.candidate
    }
    score <- sum(seq_scores[nds])
    if (score>maxScore){
      maxScore <- score
      maxClique <- nds
    }
    runtime <- Sys.time()-start.time
  }
  return(list('clique'=maxClique,'score'=maxScore))
}


clq <- max_clique_SA(g,Seqs,rc_scores,seq_scores,time.limit=1,nreps=20)
rc_table[rc_index %in% unlist(Seqs[clq$clique])]
rc_table[rc_index %in% unlist(Seqs[clq$clique]),-sum(log(P))]

profvis({SA_max_clique(g,Seqs,rc_scores,seq_scores,time.limit=1,nreps=10)})




# Annealing clique formation ----------------------------------------------





# seg_graph_to_rc_table ---------------------------------------------------




### Solving the problem for g will be a first step!

### We have more structure to our graph: the sequences represented by each node contain RC node pairs 
### Thus, we can look potentially use the graph to weight RC's and pick RC's iteratively instead of sequences
# seq_graph_to_RC_table <- function(g,VertMap,Seqs,rc_table){
#   VM <- VertMap[match(igraph::vertex_attr(g)$name,vert)]
#   VM[,degree:=igraph::degree(g)[vert]]
#   
#   rcs <- data.table('rc_index'=unique(unlist(Seqs[VM$seq])))
#   rcs[,rc_count:=c(table(unlist(Seqs[VM$seq])))[as.character(rc_index)]]
#   rcs[,rc_connectivity:=c(table(unlist(rep(Seqs[VM$seq],times=VM$degree))))[as.character(rc_index)]]
#   rcs[,rc_length:=c(table(unlist(rep(Seqs[VM$seq],times=sapply(Seqs[VM$seq],length)))))[as.character(rc_index)]/rc_count]
#        ##Average length of sequences with this rc_index
#   setkey(rcs,rc_index)
#   setkey(rc_table,rc_index)
#   rcs <- rc_table[,c('rc_index','row.node','col.node','P')][rcs]
#   setkey(rcs,rc_length,P)
#   return(rcs)
# }
# 
# 
# rcs <- seq_graph_to_RC_table(g,VertMap,Seqs,rc_table)
# 
# pick_rcix <- function(rcs){
#   rcs[rc_count==1][rc_connectivity==max(rc_connectivity)]
# }
# 
# filter_to_descendants <- function(ix,rcs){
#   ix=66903
#   RowDesc <- Row_Descendants[[toString(rcs[rc_index==ix,row.node])]] %>% unlist
#   ColDesc <- Col_Descendants[[toString(rcs[rc_index==ix,col.node])]] %>% unlist
#   rcs[row.node %in% RowDesc & col.node %in% ColDesc] %>% return
# }
# 
# 
# ### one rc, containing row.node  6089 col node 78, (6089,102) - just 4 descendants down from the microbial root 6085 - is connected to all descendants here
# ### The next most common rc has (6097:102). Col node 102 is also potentially shared with row nodes {6103,9625}
# ### with the first (6097:102) having the lowest p-value.
# 
# ### One heuristic would be:
# # 1) keep most connected rc_index - this should be a root node of every lineage - add to keeper_table
# # 2) remove our kept index and all incompatible indexes, including incompatible ancestors
# # 3) repeat until rc_table is empty
# # 4) Seqs-->max-clique on rc_indexes in keeper_table, as we may find sister clades and remove an incompatible root.
# ### With remaining 
# 
# ### Thus, we get only one lineage per sub-graph.