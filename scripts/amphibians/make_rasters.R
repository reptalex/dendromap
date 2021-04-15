# This script:
# 1. loads the ranges for the reduced set of amphibians
# 2. rasterizes their ranges while collapsing those ranges to one map per species
# 3. saves the raster to disk

library(sf)
library(raster)

amphib_ranges_reduced <- readRDS("/Users/jacobsocolar/Dropbox/Work/amphib_phylo_richness/amphib_ranges_reduced_reduced.RDS")

r <- raster(amphib_ranges_reduced, res = .1)
all_species <- unique(amphib_ranges_reduced$name_standardized)
all_species2 <- gsub(" ", "_", all_species)
for(i in 1:length(all_species)){
  print(i)
  ar <- amphib_ranges_reduced[amphib_ranges_reduced$name_standardized == all_species[i], ]
  sp_raster <- fasterize::fasterize(ar[1,], r)
  if(nrow(ar) > 1){
    for(j in 2:nrow(ar)){
      sp_raster <- sp_raster + fasterize::fasterize(ar[j,], r)
    }
    sp_raster <- sp_raster > 0
  }
  writeRaster(sp_raster, paste0("/Users/jacobsocolar/Dropbox/Work/amphib_phylo_richness/standardized_rasters/", all_species2[i], ".grd"),
              datatype="INT1S", options="COMPRESS=LZW")
}
