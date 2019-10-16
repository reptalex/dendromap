### can finding subgraphs make our life easier?
# saveRDS(G,'data/primates/big_graph')
library(data.table)
# G <- readRDS('data/primates/big_graph')
load('data/primates/primate_dendromap_guts_workspace')
sG <- igraph::clusters(G)

sG.vertices <- lapply(1:sG$no,FUN=function(a,m) names(which(m==a)),m=sG$membership)

SubGraphs <- lapply(sG.vertices,igraph::subgraph,graph=G)
lns <- sapply(sG.vertices,length)
g <- SubGraphs[[which.max(lns)]]
plot(g,vertex.size=5,vertex.label=NA,
     vertex.shape='sphere',vertex.color='black',edge.color='steelblue',
     edge.width=0.01,layout=igraph::layout_components)



SubGraphs <- SubGraphs[order(lns,decreasing = F)]
lns <- lns[order(lns,decreasing=F)]

X <- data.table('graph'=1:length(SubGraphs),
           'max_degree'=sapply(SubGraphs,FUN=function(gg) max(igraph::degree(gg))),
           'nmax'=sapply(SubGraphs,FUN=function(gg) sum(igraph::degree(gg)==max(igraph::degree(gg)))),
           'nnode'=lns)
max_possible_cliquesize <- function(g){
  degrees <- igraph::degree(g)
  tbl=table(degrees)
  x <- c(tbl)
  dd=data.table('d'=as.numeric(names(tbl)),
                'n'=rev(cumsum(rev(x))))
  dd[n>=d]
}

Cliques <- NULL
for (g in 1:length(SubGraphs)){
  base::cat(paste('\nfinding cliques in subgraph with',lns[g],'vertices'))
  clqs <- igraph::max_cliques(SubGraphs[[g]])
  Cliques <- c(Cliques,clqs)
}

plot(SubGraphs[[g]])

# Cliques <- parallel::parLapply(cl=cl,SubGraphs,igraph::max_cliques,min=2)