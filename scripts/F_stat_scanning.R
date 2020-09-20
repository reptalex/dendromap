#F-statistic scanning
rm(list=ls())
gc()
library(phylofactor)
library(dendromap)
library(parallel)
library(ggpubr)
source('R/edge_dendromap_fcns.R')
load('data/birds/bird_dendromap_workspace')

maxPval <- 0.001

rc_table <- edge_rc_table(X,row.tree,col.tree,maxPval)

# edgeMap for quick descendant calculation -------------------------------------------------------------------
rowEdgeMap <- makeEdgeMap(row.tree)
colEdgeMap <- makeEdgeMap(col.tree)


# RC_map ------------------------------------------------------------------
row.edges <- 1:Nedge(row.tree)
col.edges <- 1:Nedge(col.tree)
Row_Descendants <- rowEdgeMap[,list(list(setdiff(seq(edge,edge+n),edge))),by=edge]$V1
Col_Descendants <- colEdgeMap[,list(list(setdiff(seq(edge,edge+n),edge))),by=edge]$V1
RCmap <- edgeRCMap(rc_table,Row_Descendants,Col_Descendants)

#### simplify rc_table
rc_table <- Simplify_rc_table(rc_table,RCmap,Row_Descendants,Col_Descendants)
RCmap <- edgeRCMap(rc_table,Row_Descendants,Col_Descendants)




# scanning ----------------------------------------------------------------
rowEdgeTips <- edgeTips(row.tree)
colEdgeTips <- edgeTips(col.tree)
# pvals <- sort(unique(rc_table$P),decreasing = F)
# 
# min_pval <-  cbind(rc_table[match(RCmap$ancestor,rc_index),P],
#                    rc_table[match(RCmap$descendant,rc_index),P]) %>% apply(1,max) %>% min
# RCmap$max_P <- cbind(rc_table[match(RCmap$ancestor,rc_index),P],
#                      rc_table[match(RCmap$descendant,rc_index),P]) %>% apply(1,max)
# 
# min_ix=min(which(pvals==min_pval))
# n_pvals=100
# coarse_ix=round(seq(min_ix,length(pvals),length.out = n_pvals))
# Fstats <- numeric(n_pvals)
# 
# Lineages <- NULL
# # cl <- NULL
# t_start=Sys.time()
# for (i in coarse_ix){
#   p_thresh=pvals[i]
#   ix_thresh=rc_table[P<=p_thresh,rc_index]
#   rct <- rc_table[rc_index %in% ix_thresh]
#   rcm <- RCmap[max_P<=p_thresh]
#   
#   basal_indexes <- rcm[ancestor %in% ix_thresh & descendant %in% ix_thresh,setdiff(unique(ancestor),unique(descendant))]
#   desc_count <- rcm[ancestor %in% ix_thresh & descendant %in% ix_thresh,list(n=.N),by=ancestor]
#   basal_ixs_with_descendants <- intersect(basal_indexes,desc_count[n>1,ancestor])
#   # basal_indexes
#   # if (length(basal_ixs_with_descendants>0) & is.null(Lineages)){
#   Lineages <- lapply(basal_ixs_with_descendants,getLineage,rc_table=rct,RCmap=rcm) %>% rbindlist # method=method
#   Lineages[,lineage_size:=.N,by=lineage_id]
#   Lineages <- Lineages[lineage_size>1]
#   
#   if (length(unique(Lineages$lineage_id))>1){ ### ensure lineages don't conflict, filter accordingly
#     stats <- lineage_stats(Lineages,X,colEdgeTips,rowEdgeTips)
#     stats <- stats[order(F_stat,decreasing = T)]
#     stats <-  filter_stats(stats,Lineages,rct,rcm,Row_Descendants)
#     Lineages <- Lineages[lineage_id %in% stats$lineage_id]
#   }
#   
#   if (nrow(Lineages)>0){
#     Fstats[match(i,coarse_ix)] <- global_Fstat(Lineages,X,colEdgeTips,rowEdgeTips)
#   }
#   
#   
#   # } else {
#   #   new_indexes <- rct[P>pvals[match(i,coarse_ix)-1] & P<=p_thresh,rc_index]
#   #   new_row.edges <- rct[rc_index %in% new_indexes,row.edge]
#   #   # old_indexes <- setdiff(ix_thresh,new_indexes)
#   #   ## easy way out: re-run all basal_ix in row.tree ancestor/descendant paths of new_ix edges
#   #   row_ancs <- sapply(new_row.edges,edge_ancestors,tree=row.tree)
#   #   row_descs <- rowEdgeMap[edge %in% new_row.edges,list(seq(edge,edge+n)),by=edge]$V1 %>% unique
#   # 
#   #   affected_ix <- rct[row.edge %in% c(row_ancs,row_descs,new_row.edges),rc_index]
#   #   affected_basal_ix <- intersect(affected_ix,basal_ixs_with_descendants)
#   # 
#   #   if (length(affected_basal_ix)>0){
#   # 
#   #     if (!is.null(cl) & length(affected_basal_ix)>1){
#   #       updated_lineages <- parallel::parLapply(cl,sample(affected_basal_ix),getLineage) %>% rbindlist
#   #     } else {
#   #       updated_lineages <- lapply(affected_basal_ix,getLineage,rc_table=rct,RCmap=rcm) %>% rbindlist
#   #     }
#   #     updated_lineages[,lineage_size:=.N,by=lineage_id]
#   # 
#   #     ## only recompute stats & filter affected lineages
#   #     stats <- lineage_stats(updated_lineages,X,colEdgeTips,rowEdgeTips)
#   #     stats <- stats[order(F_stat,decreasing = T)]
#   #     stats <-  filter_stats(stats,updated_lineages,rct,rcm,Row_Descendants)
#   #     updated_lineages <- updated_lineages[lineage_id %in% stats$lineage_id]
#   # 
#   # 
#   # 
#   #     affected_lineages <- Lineages[rc_index %in% affected_ix,unique(lineage_id)]
#   # 
#   #     Lineages <- rbind(Lineages[!lineage_id %in% affected_lineages],updated_lineages)
#   #     Lineages[,lineage_size:=.N,by=lineage_id]
#   #     Lineages <- Lineages[lineage_size>1]
#   #     if (any(duplicated(Lineages[,c('row.edge','col.edge')]))){
#   #       ## this can happen if e.g. for a given row edge, col.edge 8 is a descenant of another basal row.edge with col.edge6, but
#   #       ## the same row.edge has significant col.edge 6 association and both descendants are col.edge 8, causing clean-up to assign
#   #       ## col.edge 8 to the same row.edge, producing (i,8,rc_ix=1), (i,8,rc_ix=2)
#   #       ### To resolve these disputes, we'll start off removing
#   #       # stop('duplicated stuff!')
#   # 
#   #     }
#   #     if (nrow(Lineages)>0){
#   #       if (length(unique(Lineages$lineage_id))>1){
#   #         Fstats[match(i,coarse_ix)] <- global_Fstat(Lineages,X,colEdgetips,rowEdgeTips)
#   #       } else {
#   #         Fstats[match(i,coarse_ix)] <- stats$F_stat
#   #       }
#   #     }
#   #   } else { ## else: the new_index is a singleton edge and we don't modify the Lineages table
#   #     Fstats[match(i,coarse_ix)] <- Fstats[match(i,coarse_ix)-1]
#   #   }
#   # 
#   # }
# }
# t_stop=Sys.time()
# 
# t_stop-t_start
# # parallel::stopCluster(cl)
# # rm('cl')
# 
# continent_map <- data.table('col.edge'=1:8,
#                             'continent'=c('Laurasia','Eurasia','NAmerica',
#                                           'Gondwana','Australia','SA/Africa','Africa','SouthAmerica'))
# 
# # save(list=ls(),file='data/edendromap_birds_workspace')
# 
# plot(pvals[coarse_ix],Fstats,log='x')



# more efficient way ------------------------------------------------------
stepsize=10
rowEdgeTips <- edgeTips(row.tree)
colEdgeTips <- edgeTips(col.tree)
min_pval <-  cbind(rc_table[match(RCmap$ancestor,rc_index),P],
                   rc_table[match(RCmap$descendant,rc_index),P]) %>% apply(1,max) %>% min
RCmap$max_P <- cbind(rc_table[match(RCmap$ancestor,rc_index),P],
                     rc_table[match(RCmap$descendant,rc_index),P]) %>% apply(1,max)
pvals <- sort(unique(RCmap$max_P),decreasing = F) ### need an additional trim


ps <- rc_table[,list(P=min(P),
                 rc_index=rc_index[which.min(P)]),by=row.edge]
min_P <- RCmap[ancestor %in% ps$rc_index,max_P] %>% min

pvals <- pvals[pvals>=min_P]

ix_prev <- NULL
ix <- seq(1,length(pvals),by=stepsize)
Fstats <- rep(NA,length(ix))
Lineages <- data.table()


t_start=Sys.time()
for (i in ix){
  p_thresh=pvals[i]
  ix_thresh=rc_table[P<=p_thresh,rc_index]
  new_ix <- setdiff(ix_thresh,ix_prev)
  rct <- rc_table[P<=p_thresh]
  rcm <- RCmap[max_P<=p_thresh]
  
  if (i==1){
    basal_indexes <- rcm[ancestor %in% ix_thresh & descendant %in% ix_thresh,setdiff(unique(ancestor),unique(descendant))]
    desc_count <- rcm[ancestor %in% ix_thresh & descendant %in% ix_thresh,list(n=.N),by=ancestor]
    basal_ixs_with_descendants <- intersect(basal_indexes,desc_count[n>1,ancestor])
    Lineages <- lapply(basal_ixs_with_descendants,getLins,p_thresh) %>% rbindlist
    Lineages[,lineage_size:=.N,by=lineage_id]
    Lineages <- Lineages[lineage_size>1]
  } else {
    ### here, we note that our simplify_rc_table function ensured that any new ix will:
    ### if a descendant row edge from a basal_ix, it cannot have same col.edge as prior basal_ix 
    ##   - this is because P(r,c) < P(r',c) for r' descendant of r <--> basal dominance.
    
    ### Which basal_indexes to we (re)compute?
    ### - ancestors/descendants of new_ix
    ### - note: we may even produce conflicting lineages, e.g. Rheas Gondwanan and Rheas Australian
    ### we DON'T compute: basal_ix NOT on row-tree path (root-tips) from new_ix.
    ### so... 
    ### (1) row_path_ix: Set of ix_thresh on row.tree root-tip path of new_ix
    ### (2) remove all Lineages containing row_path_ix
    ### (3) find basal_ix among row_path_ix
    ### (4) compute Lineages for basal_ix
    ### (5) rbind with old lineages
    ### (6) stats-->filter-->global F stat
    
    prev_row.edges <- rct[rc_index %in% ix_prev,unique(row.edge)]
    new_row.edges <- rct[rc_index %in% new_ix,unique(row.edge)]
    
    ## find new adges ancestral to previous edges
    new_descs <- Row_Descendants[new_row.edges]
    new_anc_edges <- sapply(new_descs,FUN=function(a,b) any (b %in% a), b=prev_row.edges) %>% new_row.edges[.]
    
    ## find previous edges ancestral to new edges
    prev_descs <- Row_Descendants[prev_row.edges]
    prev_anc_edges <- sapply(prev_descs,FUN=function(a,b) any (b %in% a), b=new_row.edges) %>% prev_row.edges[.]
    
    ## need to recompute all edges which are either equal, ancestral, or descendant to new edges
    row_path_edges <- c(intersect(prev_row.edges,new_row.edges),
                        prev_anc_edges,
                        new_anc_edges) %>% unique
    
    row_path_ix <- rct[row.edge %in% row_path_edges,unique(rc_index)]
    ## NOTE: Since clean_sisters can occasionally introduce a new row.edge not found
    ##       we'll remove all lineage_id's containing row_edges in our row_path_edges set above
    Old_Lineages <- Lineages[!row.edge %in% row_path_edges]
    
    basal_ix <- rcm[ancestor %in% row_path_ix & descendant %in% ix_thresh,setdiff(unique(ancestor),unique(descendant))]
    desc_count <- rcm[ancestor %in% row_path_ix & descendant %in% ix_thresh,list(n=.N),by=ancestor]
    basal_ixs_with_descendants <- intersect(basal_ix,desc_count[n>1,ancestor])
    Lineages <- lapply(basal_ixs_with_descendants,getLins,rc_table=rct,RCmap=rcm) %>% rbindlist
    Lineages[,lineage_size:=.N,by=lineage_id]
    Lineages <- Lineages[lineage_size>1]
    Lineages <- rbind(Old_Lineages,Lineages)
  }
  
  
  if (nrow(Lineages)>0){
    stats <- lineage_stats(Lineages,X,colEdgeTips,rowEdgeTips)
    stats <- stats[order(F_stat,decreasing = T)]
    stats <-  filter_stats(stats,Lineages,rct,rcm,Row_Descendants)
    Lineages <- Lineages[lineage_id %in% stats$lineage_id]
  }
  if (nrow(Lineages)>0){
    Fstats[match(i,ix)] <- global_Fstat(Lineages,X,colEdgeTips,rowEdgeTips)
  }
  
  ix_prev <- ix_thresh  ## we only have to recompute new_ix that affect our Lineages table
}
t_stop=Sys.time()

plot(pvals[ix],Fstats,log='x',type='o',xlab='P value cutoff',ylab='F statistic')







### this serial version may be improved... 
### it can either be brute-force parallelized with random assignment, recomputing lineages
### or we can start even coarser and iterate finer.
### this took 8.470861 mins
### Can parallelize:
### Pass rc_table, RCmap, etc. to cluster, then parLapply getLin with p_thresh input

### final output: use P_max
P_max <- pvals[ix][which.max(Fstats)]

p_thresh=P_max
ix_thresh=rc_table[P<=p_thresh,rc_index]
new_ix <- setdiff(ix_thresh,ix_prev)
rct <- rc_table[rc_index %in% ix_thresh]
rcm <- RCmap[max_P<=p_thresh]
basal_indexes <- rcm[ancestor %in% ix_thresh & descendant %in% ix_thresh,setdiff(unique(ancestor),unique(descendant))]
desc_count <- rcm[ancestor %in% ix_thresh & descendant %in% ix_thresh,list(n=.N),by=ancestor]
basal_ixs_with_descendants <- intersect(basal_indexes,desc_count[n>1,ancestor])
Lineages <- lapply(basal_ixs_with_descendants,getLins,rc_table=rct,RCmap=rcm) %>% rbindlist
Lineages[,lineage_size:=.N,by=lineage_id]
Lineages <- Lineages[lineage_size>1]
stats <- lineage_stats(Lineages,X,colEdgeTips,rowEdgeTips)
stats <- stats[order(F_stat,decreasing = T)]
stats <-  filter_stats(stats,Lineages,rct,rcm,Row_Descendants)
Lineages <- Lineages[lineage_id %in% stats$lineage_id]





# simplified: -------------------------------------------------------------
stepsize=10
ncores=7
# rowEdgeTips <- edgeTips(row.tree)
# colEdgeTips <- edgeTips(col.tree)

#### Set up P-values for scnaning
min_pval <-  cbind(rc_table[match(RCmap$ancestor,rc_index),P],
                   rc_table[match(RCmap$descendant,rc_index),P]) %>% apply(1,max) %>% min
RCmap$max_P <- cbind(rc_table[match(RCmap$ancestor,rc_index),P],
                     rc_table[match(RCmap$descendant,rc_index),P]) %>% apply(1,max)
pvals <- sort(unique(RCmap$max_P),decreasing = F) ### need an additional trim
ps <- rc_table[,list(P=min(P),
                     rc_index=rc_index[which.min(P)]),by=row.edge]
min_P <- RCmap[ancestor %in% ps$rc_index,max_P] %>% min
pvals <- pvals[pvals>=min_P]
scanned_pvals <- pvals[seq(1,length(pvals),by=stepsize)]

#### intialize cluster
cl <- parallel::makeCluster(ncores)
clusterExport(cl,varlist=setdiff(ls(),'cl'))
clusterEvalQ(cl,{library(data.table)
                 library(magrittr)})
start_time <- Sys.time()
Fscan=scanFstats(scanned_pvals,cl=cl)
stop_time <- Sys.time()
# Time difference of 1.984396 mins

### we'll refine - scanning P-values straddling the maximum

ix <- which(scanned_pvals==tst[Fstat==max(Fstat),P_thresh])
refined_pvals  <- pvals[pvals>=scanned_pvals[ix-2] & pvals<=scanned_pvals[ix+2]]

Fscan_refined <- scanFstats(refined_pvals,cl=cl)

Fscan <- rbind(Fscan,Fscan_refined)
Fscan <- Fscan[!duplicated(Fscan)]


parallel::stopCluster(cl)

ggplot(Fscan,aes(P_thresh,Fstat))+
  geom_line()+
  geom_point()+
  theme_bw()+
  scale_x_continuous(trans='log')




save(list=ls(),file='data/dendromap_Fscan_workspace')
