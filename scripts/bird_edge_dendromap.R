rm(list=ls())
gc()
library(phylofactor)
# library(dendromap)
library(parallel)
library(ggpubr)
source('Old_R/edge_dendromap_fcns.R')
load('data/birds/bird_dendromap_workspace')

set.seed(1)
start_time <- Sys.time()
dm.e <- dendromap.e(X,row.tree,col.tree,ncores=7)
stop_time <- Sys.time()
stop_time-start_time

Yhat <- predict.dendromap(dm.e)
g_obs <- plot.dendromap(dm.e,heatmap.offset = -5)
g_pred <- plot.dendromap(dm.e,Yhat,heatmap.offset = -5)
plot.new()
ggarrange(g_obs,g_pred)

ggplot(dm.e$F_scan,aes(P_thresh,Fstat))+
  geom_line()+
  geom_point()+
  theme_bw()+
  scale_x_continuous(trans='log')

Lineages <- dm.e$Lineages
rhea_id <- Lineages[row.edge==1,lineage_id]


pl <- lineage_plot(dm.e,id=rhea_id,
                   row.point.size = 7,highlight_basal=FALSE,
                   col.point.size=5,heatmap.offset = -1)
pl

save(list=ls(),file='data/edendromap_birds_workspace')
