#### Developing New algo to partition graphs based on sequence novelty scores

library(data.table)
devtools::install_github('reptalex/dendromap')
library(dendromap)
load('data/primates/primate_dendromap_guts_workspace')


# Functions ---------------------------------------------------------------
seq_score <- function(seq,scores.=scores) sum(scores[as.character(seq)])
seq_novelty_score <- function(seq,ref,scores.=scores) sum(scores[as.character(setdiff(seq,ref))])
scaffold_to_rc <- function(scaffold,Seqs.=Seqs) unique(unlist(Seqs[scaffold]))

graph_seqs <- function(g){
  igraph::vertex_attr(g)$name %>% sapply(strsplit,'_') %>%
    sapply('[',2) %>% as.numeric %>% return
}

most_novel_seq <- function(g,scaffold=NULL,Seqs.=Seqs){
  ix <- graph_seqs(g)
  ix <- setdiff(ix,scaffold)
  ### the graph vertices are indexes of our Seqs. We need to calculate novelty scores for all of those
  ### First, let's get novel rc_indexes in each
  scaffold_rc_indexes <- scaffold_to_rc(scaffold,Seqs)
  novelty_scores <- sapply(Seqs[ix],seq_novelty_score,ref=scaffold_rc_indexes)
  winner <- ix[which.max(novelty_scores)]
  # output <- list(Seqs[[winner]])
  # names(output) <- winner
  return(winner)
}

remove_incompatible_verts <- function(g,new_winner){
  vert <- paste('v_',new_winner,sep='')
  neighbors <- igraph::neighborhood(g,order = 1,nodes = vert)[[1]]$name
  return(igraph::induced_subgraph(g,neighbors))
}

check_complete <- function(g,scaffold) all(sort(paste('v_',scaffold,sep='')) ==
                                             sort(igraph::vertex_attr(g)$name))

subgraph_to_lineage <- function(sg){
  scaffold <- NULL
  while(igraph::vcount(sg)>length(scaffold)){
    new_winner <- most_novel_seq(sg,scaffold)
    scaffold <- c(scaffold,new_winner)
    sg <- remove_incompatible_verts(sg,new_winner)
  }
  return(scaffold)
}

# Scores and novelty scores -----------------------------------------------
scores <- -log(rc_table$P)
names(scores) <- rc_table$rc_index




# The algo: find most novel seq, remove incompatibles, and repeat ---------
sG <- igraph::clusters(G)
sG.vertices <- lapply(1:sG$no,FUN=function(a,m) names(which(m==a)),m=sG$membership)
SubGraphs <- lapply(sG.vertices,igraph::subgraph,graph=G)
lns <- sapply(sG.vertices,length)
g <- SubGraphs[[which.max(lns)]]
# plot(g,vertex.size=5,vertex.label=NA,
#      vertex.shape='sphere',vertex.color='black',edge.color='steelblue',
#      edge.width=0.01,layout=igraph::layout_components)

### As we iterate, we need to keep track of our scaffold for that graph.
### The scaffold will be the set of vertices, i.e. Seqs, that we fix in that graph



scaffolds <- lapply(SubGraphs,most_novel_seq)
scaffold_complete <- mapply(check_complete,SubGraphs,scaffolds) 

### pseudo-code from here:
## 0) Input: Seqs, rc_table, SubGraphs
## 1) while !all(scaffold_complete)
##     2) pick a SubGraph with incomplete scaffold
##     3) remove vertices not connected to *all* vertices in our scaffold
##     4) add remaining vertex with max_novelty_score to scaffold
##     5) repeat 3-4 until scaffold_complete=TRUE
##     6) OPTIONAL: repeat 3-5 with subgraph of removed vertices, yielding many scaffolds.
##                  Pick scaffold with largest score.
## 7) Output: One scaffold per sub-graph
## 8) convert scaffolds into rc_indexes; call list "Lineages".

## This output - Lineages of rc_indexes - can be modularly inserted in place of max_clique output.

### note: we want to keep tabs on compatibility clusters, as identifying one as "the lineage" will
### allow us to drop all incompatible alternatives within that lineage (may be faster than removal through tree traversal of descendant rc)

# making one scaffold -----------------------------------------------------

g ## our first graph

sg <- g
scaffold <- NULL
while(igraph::vcount(sg)>length(scaffold)){
  new_winner <- most_novel_seq(sg,scaffold)
  scaffold <- c(scaffold,new_winner)
  sg <- remove_incompatible_verts(sg,new_winner)
}

verts <- igraph::vertex.attributes(g)$name

cols <- rep('black',length(verts))
cols[verts %in% igraph::vertex.attributes(sg)$name] <- 'blue'




# lineages for all subgraphs ----------------------------------------------

Lineages <- lapply(SubGraphs,subgraph_to_lineage)


Lineages



# trial -------------------------------------------------------------------

ls()


dm <- dendromap(X,row.tree,col.tree,ncores=7,Pval_threshold = 0.005)
