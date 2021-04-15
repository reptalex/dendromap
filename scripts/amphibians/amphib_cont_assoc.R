# This script:
# 1. loads the reduced, harmonized amphibian maps from amphib_range_tree_harmonization.R
# 2. loads a shapefile of continents including major land-bridge islands
# 3. Extracts presence/absence of each amphib/continent combo
# 4. Saves the presence/absence matrix (RDS save of dataframe object)
library(sf)

amphib_ranges_reduced <- readRDS("/Users/jacobsocolar/Dropbox/Work/amphib_phylo_richness/amphib_ranges_reduced.RDS")

setwd("/Users/JacobSocolar/Dropbox/Work/Diversity_accum")
load("continents.Rdata")
geometry <- lapply(continents, function(x) st_as_sf(x))
conts <- do.call(rbind, geometry)
conts$id <- names(continents)
conts2 <- st_transform(conts, st_crs(amphib_ranges_reduced))

amphib_conts <- st_intersects(amphib_ranges_reduced, conts2, sparse = F)

amphib_conts <- as.data.frame(amphib_conts)
names(amphib_conts) <- names(continents)
amphib_conts <- amphib_conts[,2:6]
amphib_conts$species <- amphib_ranges_reduced$name_standardized
sum(duplicated(amphib_conts$species))

# collapse over the duplicated rows
amphib_conts2 <- amphib_conts
for(i in 1:nrow(amphib_conts)){
  ac2 <- amphib_conts[amphib_conts$species == amphib_conts$species[i],]
  ac3 <- amphib_conts2[amphib_conts2$species == amphib_conts$species[i],]
  if(nrow(ac2) > 1 & nrow(ac3) > 1){
    amphib_conts2 <- amphib_conts2[-which(amphib_conts2$species == ac2$species[1]), ]
    newrow <- c(colSums(ac2[1:5]) > 0, ac2$species[1])
    amphib_conts2 <- rbind(amphib_conts2, newrow)
  }
}

# convert to numeric
for(i in 1:5){
  amphib_conts2[,i] <- as.numeric(as.logical(amphib_conts2[,i]))
}
head(amphib_conts2)

# Ensure all rowSums > 0
sum(rowSums(amphib_conts2[,1:5])>0)
dim(amphib_conts2)

saveRDS(amphib_conts2, file = "/Users/jacobsocolar/Dropbox/Work/amphib_phylo_richness/amphib_conts_standardized.RDS")
