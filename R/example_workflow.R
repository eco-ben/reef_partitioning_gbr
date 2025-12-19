# Install package using GitHub repo and specific template branch
remotes::install_github("VHallerBull/ReefPartitionAllenAtlas", ref="package-template")

library(ReefPartitionAllenAtlas)

library(sf)
library(terra)
library(h3)
library(tidyverse)
library(igraph)

set.seed(123)

ID <- "23055101104" # One Tree Island Reef

# Load reef polygons data
reef_polygons <- st_read("../canonical-reefs/output/rrap_canonical_2025-01-10-T13-28-06.gpkg")

# Load bathymetry and geomorphic habitat data
# ! Ensure that these two raster datasets use the same CRS !
bathy_raster <- rast("../../reef_partitioning/Data/required_data/bathy/Mackay-Capricorn/MC_wg84.tif")
habitat_raster <- rast("../../reef_partitioning/Data/required_data/GBR10 GBRMP Geomorphic.tif")

# Must ensure all input data are in the same CRS (instead of GDA2020)
reef_polygons <- st_make_valid(st_transform(reef_polygons, crs = st_crs(habitat_raster)))

habitat_categories <- c(15, 22, 14, 21) # reef crest, reef slope, outer reef flat, sheltered reef slope

# Define hexagon resolution and unit for outputting hexagon area (same as function defaults)
hex_resolution <- 12
unit <- "km2"

# Extract pixel data for target reef - OTI reef
# This will select pixels from `habitat_raster` that cover `target_reef` and
# extract the depth values that cover these pixels as well.
target_reef <- reef_polygons[reef_polygons$UNIQUE_ID == ID, ]
pixel_data <- extract_pixel_points(
    target_reef,
    habitat_raster,
    bathy_raster,
    habitat_categories,
    hex_resolution = hex_resolution,
    unit = unit,
    additional_variable_name = "depth",
    output_epsg = 3112,
    resample_method = "bilinear"
)
pixel_data$UNIQUE_ID <- target_reef$UNIQUE_ID

# Cluster the pixels into site clusters using constr.clust algorithm with a
# minimum spanning tree input using default settings
mst_hclust_pixels <- cluster_reef_pixels(
    pixel_data,
    habitat_clustering_function = function(x) {
        mst <- prepare_mst(x)
        clusters <- constrained_hclust(x, as_edgelist(mst))
        clusters
    }
)
# Collate clustered site pixels into site polygons
mst_hclust_sites <- clustered_pixels_to_polygons(mst_hclust_pixels)

# Cluster the pixels into site clusters using skater algorithm with minimum spanning
# tree (mst is created within reef_skater())
skater_pixels <- cluster_reef_pixels(
    pixel_data,
    habitat_clustering_function = reef_skater,
    parallelisation = "Windows"
)
skater_sites <- clustered_pixels_to_polygons(skater_pixels)