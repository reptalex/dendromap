rm(list=ls())
gc()
source('R/edge_dendromap_fcns.R')
load('data/birds/bird_dendromap_workspace')
m <- nrow(X)
n <- ncol(X)

set.seed(1)

# edge rc_table -----------------------------------------------------------
rc_table <- edge_rc_table(X,row.tree,col.tree)

# edgeMap for quick descendant calculation -------------------------------------------------------------------
rowEdgeMap <- makeEdgeMap(row.tree)
colEdgeMap <- makeEdgeMap(col.tree)

# Rhea case study ----------------------------------------------------
############ CASE STUDY: RHEAS ###############
## Rheas: edge 1 (Rhea) and 106 (non-Rhea) Need: 10
rc_table[row.edge==1]
# row.edge col.edge     stat   rank            P rc_index
# 1:        1        4 2.016367 3409.5 7.455945e-06     1435 #Gondwana
# 2:        1        6 2.239922 2681.5 5.705469e-06     1884 #SA/Af
# 3:        1        8 4.217850  615.5 8.964294e-07     3268 #SA

### SA/Africa is more significant because only Emus/Cassowaries live on Australia. 
### We'd pick Gondwana if and only if we ID BOTH an Africa AND 

rheas <- row.tree$tip.label[phangorn::Descendants(row.tree,8106,'tips')[[1]]]
rhea.tree <- drop.tip(row.tree,setdiff(row.tree$tip.label,rheas))

par(mfrow=c(1,2))
plot(rhea.tree)
edgelabels(1+1:(ape::Nedge(rhea.tree)))
plot(col.tree)
edgelabels()

rc_table[row.edge %in% rowEdgeMap[edge==1,seq(edge,edge+n)]]
rc_table[row.edge==1]

rc_table[row.edge==1+95] ## cassowaries & Emus - Australia
rc_table[row.edge==1+1]  ## tinamous - south america
rc_table[row.edge==1+94]  ## tinamous - south america

focal_edge=1
descendants <- rowEdgeMap[edge==focal_edge,setdiff(seq(edge,edge+n),edge)]
rhea_table <- rc_table[row.edge %in% c(descendants,focal_edge)]


# RC_map ------------------------------------------------------------------
row.edges <- 1:length(phylofactor::getPhyloGroups(row.tree))
col.edges <- 1:length(phylofactor::getPhyloGroups(col.tree))
Row_Descendants <- rowEdgeMap[,list(list(setdiff(seq(edge,edge+n),edge))),by=edge]$V1
Col_Descendants <- colEdgeMap[,list(list(setdiff(seq(edge,edge+n),edge))),by=edge]$V1
RCmap <- edgeRCMap(rc_table,Row_Descendants,Col_Descendants)

# Rhea analysis -----------------------------------------------------------


# ancs <- unique(rhea_map$ancestor)
# desc <- unique(rhea_map$descendant)
# basal_edges <- setdiff(ancs,desc)

focal_edge=1
descendants <- rowEdgeMap[edge==focal_edge,setdiff(seq(edge,edge+n),edge)]
rhea_table <- rc_table[row.edge %in% c(descendants,focal_edge)]
rhea_map <- edgeRCMap(rhea_table,Row_Descendants,Col_Descendants)
### old route: joinability graph partitioning
# rhea_nds <- unique(rhea_map[terminal==TRUE,descendant])
# rhea_sqs <- lapply(rhea_nds,rc_seqs,rhea_map) %>% unlist(recursive=FALSE) %>% unique
# rhea_joinables <- find_joinables(rhea_sqs,rhea_table,Row_Descendants,Col_Descendants)
rhea_lineage <- greedyColTree(rhea_table,rhea_map,rc_table,row.tree)

# > rhea_lineage
# row.edge col.edge     stat   rank            P rc_index lineage_id
# 1:        1        8 4.217850  615.5 8.727693e-07     3385       3385
# 2:       96        5 1.992751 3539.0 6.295774e-06     1450       3385
# 3:        2        6 2.508307 2107.0 3.638959e-06     2375       3385
### looks good - Gondwanan Rheas, South American tinamous, Australian Emus/Cassowaries
### We want to make sure this table is found in the bird dendromap

# choosing edges to input to greedyColTree -----------------------------------------
ancs <- unique(RCmap$ancestor)
descs <- unique(RCmap$descendant)
basal_indexes <- setdiff(ancs,descs)

desc_count <- RCmap[,list(n=.N),by=ancestor]

basal_ixs_with_descendants <- intersect(basal_indexes,desc_count[n>1,ancestor])

# lineages <- lapply(basal_indexes,getLineage,
#                         RCmap,rc_table,row.tree,Row_Descendants)

cl <- parallel::makeCluster(7)

parallel::clusterExport(cl,varlist=c('incompatible_descendants','clean_sisters','greedyColTree','getLineage','edge_ancestors',
                                     'RCmap','rc_table','Row_Descendants','Col_Descendants','row.tree','col.tree'))
parallel::clusterEvalQ(cl,library(data.table))
parallel::clusterEvalQ(cl,library(magrittr))
Lineages <- parallel::parLapply(cl,sample(basal_ixs_with_descendants),getLineage) %>% rbindlist
parallel::stopCluster(cl)
rm('cl')

Lineages[,n:=.N,by=lineage_id]
Lineages <- Lineages[n>1]



# Filtering lineages ------------------------------------------------------

colEdgeTips <- edgeTips(col.tree)
rowEdgeTips <- edgeTips(row.tree)

### get fit stats
lineage_ids <- Lineages[,list(n=.N),by=lineage_id][n>1]$lineage_id
stats <- lineage_stats(Lineages,X,colEdgeTips,rowEdgeTips)
stats <- stats[order(F_stat,decreasing = T)]

filtered_stats <- filter_stats(stats)

## check for rheatites - 1435 - and suboscines - 3632
c(1435,3632) %in% filtered_stats$lineage_id


Lineages <- Lineages[lineage_id %in% filtered_stats$lineage_id]


setkey(Lineages,lineage_id,row.edge)
setkey(filtered_stats,lineage_id)

Lineages <- filtered_stats[,c('lineage_id','F_stat')][Lineages]
setkey(Lineages,row.edge)
Lineages[,row.node:=edge2node(row.edge,row.tree)]
Lineages[,col.node:=edge2node(col.edge,col.tree)]
Lineages <- Lineages[order(F_stat,decreasing = T)]
Lineages[,Lineage:=match(lineage_id,unique(lineage_id))]
object <- list('Lineages'=Lineages,
               'Data'=X,
               'row.tree'=row.tree,
               'col.tree'=col.tree,
               'et.c'=colEdgeTips,
               'et.r'=rowEdgeTips,
               'desc.r'=Row_Descendants,
               'desc.c'=Col_Descendants)
class(object) <- 'dendromap'



plot.dendromap(object)

Yhat <- predict.dendromap(object)

g_obs <- plot.dendromap(object)
g_pred <- plot.dendromap(object,y=Yhat)

ggpubr::ggarrange(g_obs,g_pred)
