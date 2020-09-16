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

# demo --------------------------------------------------------------------

load('data/edendromap_birds_workspace')
stats  # Lineage F-statistics, sorted from highest to lowest

# i=55 gets the Rheas
i=55
lineage = Lineages[lineage_id==stats$lineage_id[i]]

#    row.edge col.edge     stat rank            P rc_index lineage_id lineage_size
# 1:    13270        6 2.419963 2294 4.473644e-06     7883       6689            2
# 2:    13274        7 4.993270  398 6.057702e-07     9156       6689            2

spp <- spp_lineage(lineage)
spp
# $Gondwana
# [1] "Nothura"      "Taoniscus"    "Nothoprocta"  "Rhynchotus"   "Eudromia"     "Tinamotis"    "Nothocercus" 
# [8] "Crypturellus" "Tinamus"      "Struthio"     "Dromaius"     "Casuarius"    "Rhea"        
# 
# $Australia
# [1] "Dromaius"  "Casuarius" "Rhea"     
# 
# $SouthAmerica
# [1] "Nothura"      "Taoniscus"    "Nothoprocta"  "Rhynchotus"   "Eudromia"     "Tinamotis"    "Nothocercus" 
# [8] "Crypturellus" "Tinamus"    
