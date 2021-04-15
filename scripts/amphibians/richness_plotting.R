# Plot the species richness of a list of amphibians
library(raster)

# Get an example species list:
amphib_ranges_reduced <- readRDS("/Users/jacobsocolar/Dropbox/Work/amphib_phylo_richness/amphib_ranges_reduced.RDS")
species_list <- unique(amphib_ranges_reduced$name_standardized[1:50])

# Pointer to the locations of the single-species rasters from make_rasters.R
raster_location <- "/Users/jacobsocolar/Dropbox/Work/amphib_phylo_richness/standardized_rasters/"

# Define the get_richness_raster function
source("/Users/jacobsocolar/Dropbox/Work/Code/amphibian_dendromap/get_richness_raster.R")

# get an example richness raster
ex_richness <- get_richness_raster(species_list, raster_location)

# plot it
raster::plot(ex_richness)

##### make it look nicer: ######
# load up the continents as a raster
conts_raster <- raster::raster("/Users/jacobsocolar/Dropbox/Work/amphib_phylo_richness/continents_raster.grd")
ex_richness[is.na(conts_raster)] <- NA
# Now we can distinguish continent outlies from oceans
plot(ex_richness)

# reproject to a better (? YMMV) projection
crs <- "+proj=robin +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m"
robinson_richness <- projectRaster(ex_richness, crs=crs)
plot(robinson_richness)

# mask out duplicated land areas
library(geosphere)
e <- as(extent(ex_richness), "SpatialPolygons")
crs(e) <- crs(ex_richness)
e <- makePoly(e)
re <- spTransform(e, crs)

robinson_richness_masked <- mask(robinson_richness, re)
raster::plot(robinson_richness_masked, axes=F, box = F,
             colNA = "white", col = viridis::viridis(100))

# soften the edges a bit
rrm_smooth <- raster::focal(robinson_richness_masked, w=matrix(1, 9, 9))/81

raster::plot(rrm_smooth, axes=F, box = F, col = viridis::viridis(100))

# try a black background
raster::plot(rrm_smooth, axes=F, box = F,
             colNA = "black", col = viridis::viridis(100))

# overlay white so that it looks like a robinson-projection map
edgemask <- fasterize::fasterize(sf::st_as_sf(re), rrm_smooth, background = 2)
edgemask[edgemask==1] <- NA
plot(edgemask, col = "white", add = T, legend = F)
