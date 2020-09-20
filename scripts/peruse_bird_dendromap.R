library(ape)
library(data.table)
library(magrittr)


# functions ---------------------------------------------------------------

get_spp <- function(x,row.tree.=row.tree,rowEdgeTips.=rowEdgeTips) row.tree$tip.label[rowEdgeTips[edge==x,seq(min,max)]]
spp_lineage <- function(lineage,continent_map.=continent_map,genera=TRUE){
  spp <- lapply(lineage$row.edge,get_spp)
  if (genera){
    spp2gen <- function(spp) sapply(spp,strsplit,'_') %>% sapply(getElement,1) %>% unique
    spp <- lapply(spp,spp2gen)
  }
  names(spp) <- continent_map[match(lineage$col.edge,col.edge),continent]
  return(spp)
}


continent_map <- data.table('col.edge'=1:8,
                            'continent'=c('Laurasia','Eurasia','NAmerica',
                                          'Gondwana','Australia','SA/Africa','Africa','SouthAmerica'))

# demo --------------------------------------------------------------------

load('data/edendromap_birds_workspace')
stats <- Lineages[,list(F_stat=unique(F_stat)),by=lineage_id][order(F_stat,decreasing=T)]
stats  # Lineage F-statistics, sorted from highest to lowest
rowEdgeTips <- edgeTips(row.tree)


i=1
i=i+1
lineage = Lineages[lineage_id==stats$lineage_id[i]]
# lineage
lid <- lineage$lineage_id[1]
lid
lineage_plot(dm.e,id=lid,
             row.point.size = 7,highlight_basal=FALSE,
             col.point.size=5,heatmap.offset = -1)

spp <- spp_lineage(lineage)
spp
