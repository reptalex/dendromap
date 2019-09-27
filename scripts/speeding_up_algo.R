library(dendromap)
library(tictoc)
library(lineprof)
library(profvis)
library(parallel)
set.seed(1)
m=1000
n=10
row.tree <- rtree(m)
col.tree <- rtree(n)
S <- treeSim(5,row.tree,col.tree,col.node=n+1)
eta <- S$W %*% (10*S$D) %*% t(S$V)
X <- eta+matrix(rnorm(m*n),nrow=m)
clrinv <- function(x) exp(x)/sum(exp(x))
rmlt <- function(p,lambda=5e3) rmultinom(1,rpois(1,lambda),prob = p)
N <- apply(X,2,clrinv) %>% apply(2,rmlt)
rownames(N) <- row.tree$tip.label
colnames(N) <- col.tree$tip.label

base::cat(paste('\nMaking Nodemaps'))
row.nodemap <- dendromap:::makeNodeMap(row.tree)
col.nodemap <- dendromap:::makeNodeMap(col.tree)

# base::cat(paste('\nMaking RC table with',n_sim,'null simulations'))
rc_table <- makeRCtable(X,row.tree,col.tree)
rc_table <- rc_table[P<max(P)]
rc_table[,fdr:=p.adjust(P,'fdr')]
Pval_threshold <- 0.05
rc_table <- rc_table[P<=Pval_threshold]
RCmap <- makeRCMap(rc_table,row.nodemap,col.nodemap)

# tic()
# Lineages <- find_lineages(RCmap,rc_table,row.nodemap,col.nodemap)
# toc()
# ## 30s
# 
# cl <- makeCluster(3)
# clusterEvalQ(cl,library(dendromap))
# clusterExport(cl,setdiff(ls(),'cl'))
# clusterExport(cl,'find_lineages2')
# tic()
# LL=find_lineages2(RCmap,rc_table,row.nodemap,col.nodemap,cl)
# toc()
# ## 12s
# stopCluster(cl)
# rm('cl')


# RCmap -------------------------------------------------------------------

### two important things here:
# 1) reduce time, though that's less urgent than:
# 2) reduce size of RCmap - minimize to only nearest-ancestor

row.nodes <- unique(rc_table$row.node)
col.nodes <- unique(rc_table$col.node)
Row_Descendants <- lapply(row.nodes,getIndexSets,row.nodemap) %>%
  lapply(FUN=function(x,a) lapply(x,intersect,a),a=row.nodes)
Col_Descendants <- lapply(col.nodes,getIndexSets,col.nodemap) %>%
  lapply(FUN=function(x,a) lapply(x,intersect,a),a=col.nodes)

names(Row_Descendants) <- row.nodes
names(Col_Descendants) <- col.nodes

terminalRowNodes <- row.nodes[sapply(Row_Descendants,FUN=function(x) length(unlist(x)))==0]
terminalColNodes <- row.nodes[sapply(Row_Descendants,FUN=function(x) length(unlist(x)))==0]
## we don't have to look for descendants of these row/col nodes

ix <- which((!rc_table$row.node %in% terminalRowNodes) &
              (!rc_table$col.node %in% terminalColNodes))
maps <- lapply(ix,makeDescendantTable,
               rc_table,
               terminalRowNodes,
               terminalColNodes,
               Row_Descendants,
               Col_Descendants)
RCmap <- rbindlist(maps)
RCmap[,terminal:=!(descendant %in% ancestor)]
setkey(RCmap,descendant,ancestor)

j <- table(RCmap[,descendant]) %>% which.max %>% names %>% as.numeric

RCmap[descendant==j]


# rc_seqs -------------------------------------------------

### even faster will be an rc_seq algo that, if find 1-2-3, removes 1-3 from table.
# used a new rc_seqs that only uses the terminal nodes, NOT their rows:

nds <- unique(RCmap[terminal==TRUE,descendant])

ix <- which(RCmap$terminal)
tic()
Seqs <- lapply(ix,rc_seqs,RCmap) %>% unlist(recursive=FALSE)
toc()
## 30s

cl <- makeCluster(3)
clusterEvalQ(cl,library(dendromap))
clusterExport(cl,setdiff(ls(),'cl'))
nds <- unique(RCmap[terminal==TRUE]$descendant)
tic()
Seqs2 <- parLapply(cl,nds,rc_seqs,RCmap) %>% unlist(recursive=FALSE)
toc()
## 7.8s --- WAAY FASTER! There are ~1000 nodes vs. 3000 indexes. Plus parallelization = this is good!

sqs <- sapply(Seqs2,FUN=function(a,b) b %in% a, b=j) %>% which
Seqs2[sqs]

Seqs <- Seqs2

# check_joinable ----------------------------------------------------------
n <- length(Seqs)
tbl <- data.table('seq1'=rep(1:(n-1),times=(n-1):1),key='seq1')
tbl[,seq2:=(seq1+1):n,by=seq1]  ## this table has all pairwise indexes we need to check.
### We might be able to better index these, perhaps as far back as Seqs

############# ORIGINAL 
tic()
joinability <- parApply(cl=cl,t(tbl),2,FUN=function(x,s,rc,r,c) check_joinable(x[1],x[2],s,rc,r,c),
                        s=Seqs,rc=rc_table,r=row.nodemap,c=col.nodemap)
toc()
## 9.4s

stopCluster(cl)
rm('cl')

tic()
tic()
joinability <- apply(t(tbl),2,FUN=function(x,s,rc,r,c) check_joinable(x[1],x[2],s,rc,r,c),
                     s=Seqs,rc=rc_table,r=row.nodemap,c=col.nodemap)
toc()
## 23.7s

#### Now to see if we can speed up check_joinable
profvis({apply(t(tbl[1:1000]),2,FUN=function(x,s,rc,r,c) check_joinable(x[1],x[2],s,rc,r,c),
               s=Seqs,rc=rc_table,r=row.nodemap,c=col.nodemap)})

ixs <- unique(unlist(Seqs))

tic()
joinable <- apply(t(tbl),2,FUN=function(x,s,rc,r,c) cj(x[1],x[2],s,rc,r,c),
      s=Seqs,rc=rc_table,r=Row_Descendants,c=Col_Descendants)
toc()
#6s - SWEEET

profvis({apply(t(tbl),2,FUN=function(x,s,rc,r,c) cj(x[1],x[2],s,rc,r,c),
               s=Seqs,rc=rc_table,r=Row_Descendants,c=Col_Descendants)})

length(unique(rc_table[rc_index %in% ixs,row.node])) #132
length(unique(rc_table[rc_index %in% ixs,col.node])) #9
nrow(tbl)                                            #8256...
## We only need to compute descendants for row tree 132 times, and for col tree 9 times
## instead, the above algo has us compute it about 16,500 times!!!

## can change check_joinable to take as input 


### check-joinability is currently a bottleneck - we need to speed this up bigtime.
### One possibility: dynamically construct a reference for seqs using only their roots
  ## - for seqs with the same root, we calculate descendants only once, as getIndexSets takes the most time
### to subset tbl into seqs with possible compatibility.

### Must make sure our new output is the same as the old output


############# NEW
## 

# Cliques ----------------------------------------------------------------

tbl[,joinability:=joinability]

seq.indexes <- unique(unlist(tbl[,c('seq1','seq2')]))
nverts <- length(seq.indexes)
VertMap <- data.table('seq'=seq.indexes,'vert'=paste('v',seq.indexes,sep='_'),key='vert')

joinables <- tbl[joinability==TRUE]
nedge <- nrow(joinables)

joinable.edges <- split(joinables[,c('seq1','seq2')],seq(nrow(joinables))) %>% unlist
joinable.edges <- VertMap[match(joinable.edges,seq),vert]

########## ORIGINAL
G <- igraph::make_graph(joinable.edges,directed = F)  ##graph 
tic()
  joinable.seqs <- igraph::max_cliques(G,min=2)
  joinable.seqs <- joinable.seqs[order(sapply(joinable.seqs,length),decreasing = F)]
toc()


######### NEW 
clqs <- function(nd,G.=G,sG.=sG){
  nds <- names(which(sG$membership==nd))
  igraph::subgraph(G,nds) %>% igraph::max_cliques(min=2) %>%
    return()
}
tic()
sG <- igraph::clusters(G)            

jsqs=lapply(1:sG$no,clqs)
### we can parallelize maxClique over subgraphs

########## FIND MAXIMAL CLIQUES
### below can be parallelized by defining subsets (throw out disconnected nodes)
# joinable.seqs <- igraph::cliques(G,min=2) ### now we need to remove elements that are strict subsets of other sets
joinable.seqs <- igraph::max_cliques(G,min=2)
joinable.seqs <- joinable.seqs[order(sapply(joinable.seqs,length),decreasing = F)]


lns <- sapply(rev(joinable.seqs),length)
if (length(joinable.seqs)>10){
  big.lns <- paste(lns[1:5],collapse=',')
  small.lns <- paste(lns[(length(joinable.seqs)-4):length(joinable.seqs)],collapse=',')
} else {
  small.lns <- paste(lns,collapse=',')
  big.lns <- NULL
}
base::cat(paste('\nFound',length(joinable.seqs),'joinable cliques of sizes',big.lns,'...',small.lns))
base::cat('\nRemoving subsets to find maximal cliques')

### below can be parallelized, but now is unncessary b.c. we're using max_cliques?
# for (i in 1:(length(joinable.seqs)-1)){
#   check.subset <- any(sapply(joinable.seqs[(i+1):length(joinable.seqs)],
#                              FUN=function(a,b) all(b %in% a),b=joinable.seqs[[i]]))
#   if (any(check.subset)){
#     joinable.seqs[[i]] <- NA
#   }
# }
found.cliques <- !sapply(joinable.seqs,FUN=function(x) all(is.na(x)))
base::cat(paste('\nThere are',sum(found.cliques),'cliques remaining'))
joinable.seqs <- joinable.seqs[found.cliques]

## add disconnected RC's
all.RC.cliques <- c(as.list(setdiff(VertMap$vert,joinable.edges)),
                    lapply(joinable.seqs,names))

all.RC.cliques <- lapply(all.RC.cliques,FUN=function(v,VertMap) VertMap[match(v,vert),seq],
                         VertMap)

get_rc_index <- function(x,Seqs.=Seqs) unique(unlist(Seqs[x]))
RCSeqs <- lapply(all.RC.cliques,get_rc_index)
