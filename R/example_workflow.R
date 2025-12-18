# These packages will likely go into our package's dependencies when the collated package is working
library(sf)
library(terra)
library(h3)
library(sfnetworks)
library(igraph)
library(spdep)
library(adespatial)
library(tidyverse)
library(exactextractr)
library(furrr)

plan(multisession, workers = 6)

# library(ReefPartition) # after package has been built
# As package is not ready for building we need to source each function
source("./R/pixel_extraction.R")
source("./R/cluster_reef_pixels.R")
source("./R/clustering_algorithms/constrained_hclust.R")
source("./R/clustering_algorithms/prepare_cluster_edges.R")
source("./R/clustering_algorithms/reef_skater.R")
source("./R/site_postprocessing.R")

# Set seed for reproducibility
set.seed(123)

ID <- "23055101104" # One Tree Island Reef

reef_polygons <- st_read("../canonical-reefs/output/rrap_canonical_2025-01-10-T13-28-06.gpkg")

bathy_raster <- rast("../../reef_partitioning/Data/required_data/bathy/Mackay-Capricorn/MC_wg84.tif")
habitat_raster <- rast("../../reef_partitioning/Data/required_data/GBR10 GBRMP Geomorphic.tif")

# Must ensure all input data are in the same CRS
reef_polygons <- st_transform(reef_polygons, crs = st_crs(habitat_raster)) %>%
    st_make_valid()

habitat_categories <- c(15, 22, 14, 21) # reef crest, reef slope, outer reef flat, sheltered reef slope

# Define hexagon resolution and unit for outputting hexagon area (same as function defaults)
hex_resolution <- 12
unit <- "km2"

# Extract pixel data for target reef - OTI reef
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

pixel_data$hexagons <- h3_to_geo_boundary_sf(pixel_data$h3_index)
pixel_data$depth_extract <- abs(exact_extract(bathy_raster, pixel_data$hexagons, "mean"))
pixel_data <- select(pixel_data, !hexagons)

# Test approach with MST constr.hclust
mst_hclust_pixels <- cluster_reef_pixels(
    pixel_data,
    habitat_clustering_function = function(x) {
        mst <- prepare_mst(x)
        constrained_hclust(x, as_edgelist(mst))
    }
)
mst_hclust_sites <- clustered_pixels_to_polygons(mst_hclust_pixels)

# Test approach with triangulated constr.hclust
tri_hclust_pixels <- cluster_reef_pixels(
    pixel_data,
    habitat_clustering_function = function(x) {
        tri_edges <- prepare_tri_edges(x)
        constrained_hclust(x, tri_edges)
    }
)
tri_hclust_sites <- clustered_pixels_to_polygons(tri_hclust_pixels)

# Test approach with high weight constr.hclust
high_weight_hclust_pixels <- cluster_reef_pixels(
    pixel_data,
    habitat_clustering_function = function(x) {
        tri_edges <- prepare_tri_edges(x)
        constrained_hclust(x, tri_edges, alpha = 0.7)
    }
)
high_weight_hclust_sites <- clustered_pixels_to_polygons(high_weight_hclust_pixels)

# Test approach with skater (as in previous functions)
skater_pixels <- cluster_reef_pixels(
    pixel_data,
    habitat_clustering_function = function(x) reef_skater(x)
)
saveRDS(skater_pixels, "OTIReef_skater.rds")
skater_pixels <- readRDS("OTIReef_skater.rds")
skater_sites <- clustered_pixels_to_polygons(skater_pixels)

extract_habitat <- function(id) {
    one <- str_extract(as.character(id), "_\\d{2}_")
    two <- str_extract(one, "\\d{2}")
}

skater_compactness <- skater_sites %>%
    st_transform(crs = 3112) %>%
    mutate(habitat = extract_habitat(site_id)) %>%
    mutate(conv_hull = st_convex_hull(geometry), minbound_circle = st_minimum_bounding_circle(geometry)) %>%
    mutate(conv_hull_area = st_area(conv_hull), circle_area = st_area(minbound_circle), site_area = st_area(geometry)) %>%
    mutate(prop_conv_hull = site_area / conv_hull_area, prop_circle = site_area / circle_area)

skater_pixels$hexagons <- h3_to_geo_boundary_sf(skater_pixels$h3_index)
skater_pixels$depth <- abs(exact_extract(bathy_raster, skater_pixels$hexagons, "mean"))

skater_variance <- skater_pixels %>%
    group_by(habitat, site_id) %>%
    summarise(depth_sd = sd(depth), .groups = "keep") %>%
    group_by(habitat) %>%
    summarise(mean_habitat_sd = mean(depth_sd), variance_habitat_sd = var(depth_sd))

mean_skater_slope_sd <- as.numeric(mean(skater_variance[skater_variance$habitat == 22, ]$mean_habitat_sd))
mean_skater_convhull_compact <- as.numeric(mean(skater_compactness[skater_compactness$habitat == 22, ]$prop_conv_hull, na.rm=TRUE))
mean_skater_circle_compact <- as.numeric(mean(skater_compactness[skater_compactness$habitat == 22, ]$prop_circle, na.rm = TRUE))

compare_skater_slopes <- function(param) {
    alpha <- param$alpha
    beta <- param$beta
    method <- param$method

    hclust_pixels <- cluster_reef_pixels(
        pixel_data,
        habitat_clustering_function = function(x) {
            mst <- prepare_mst(x)
            constrained_hclust(x, as_edgelist(mst), alpha = alpha, beta = beta, method = "average")
        }
    )
    hclust_sites <- clustered_pixels_to_polygons(hclust_pixels)

    hclust_compactness <- hclust_sites %>%
        st_transform(crs = 3112) %>%
        mutate(habitat = extract_habitat(site_id)) %>%
        mutate(conv_hull = st_convex_hull(geometry), minbound_circle = st_minimum_bounding_circle(geometry)) %>%
        mutate(conv_hull_area = st_area(conv_hull), circle_area = st_area(minbound_circle), site_area = st_area(geometry)) %>%
        mutate(prop_conv_hull = site_area / conv_hull_area, prop_circle = site_area / circle_area)

    hclust_variance <- hclust_pixels %>%
        group_by(habitat, site_id) %>%
        summarise(depth_sd = sd(depth), .groups = "keep") %>%
        group_by(habitat) %>%
        summarise(mean_habitat_sd = mean(depth_sd), variance_habitat_sd = var(depth_sd))

    hclust_depth_sd <-  as.numeric(mean(hclust_variance[hclust_variance$habitat == 22, ]$mean_habitat_sd))
    hclust_convhull_compact <-  as.numeric(mean(hclust_compactness[hclust_compactness$habitat == 22, ]$prop_conv_hull))
    hclust_circle_compact <- as.numeric(mean(hclust_compactness[hclust_compactness$habitat == 22, ]$prop_circle))

    # skater_depth_sd_diff <- mean_skater_slope_sd - hclust_depth_sd
    # skater_convhull_diff <- mean_skater_convhull_compact - hclust_convhull_compact
    # skater_circle_diff <- mean_skater_circle_compact - hclust_circle_compact
    
    results <- data.frame(
        hclust_depth_sd = hclust_depth_sd,
        hclust_convhull_compact = hclust_convhull_compact,
        hclust_circle_compact = hclust_circle_compact
    )

    return(list(results = results, clusters = hclust_pixels))
}

alphas <- seq(0, 1, 0.01)
betas <- -1

params <- expand.grid(alpha = alphas, beta = betas)
params$perm <- 1:nrow(params)

params_list <- split(params, params$perm)

param_results <- future_map(
    params_list,
    compare_skater_slopes,
    .progress = TRUE
)

param_clusters <- lapply(param_results, function(x) x$clusters)
param_clusters <- lapply(seq_along(param_clusters), function(x_ind) {
    x <- param_clusters[[x_ind]]
    x[, "alpha"] <- params_list[[x_ind]]$alpha
    x[, "beta"] <- params_list[[x_ind]]$beta
    x
})
param_clusters <- do.call(rbind,param_clusters)

param_values <- lapply(param_results, function(x) x$results)
param_values <- do.call(rbind, param_values)

params <- cbind(params, param_values)
params$alpha <- as.character(params$alpha)
params$beta <- as.character(params$beta)

# Comparing results to One Tree Original Skater results
old_skater <- st_read("../../reef_partitioning/Data/OneTree62500_2025_07_28/OneTree62500_Polygon_Geometry.gpkg")
old_skater_compact <- st_transform(old_skater, crs = 3112) %>%
    mutate(conv_hull = st_convex_hull(geom), minbound_circle = st_minimum_bounding_circle(geom)) %>%
    mutate(conv_hull_area = st_area(conv_hull), circle_area = st_area(minbound_circle), site_area = st_area(geom)) %>%
    mutate(prop_conv_hull = as.numeric(site_area / conv_hull_area), prop_circle = as.numeric(site_area / circle_area))

old_skater_convhull_slope <- mean(
    old_skater_compact[old_skater_compact$habitat == "Reef Slope", ]$prop_conv_hull
)
old_skater_circle_slope <- mean(
    old_skater_compact[old_skater_compact$habitat == "Reef Slope", ]$prop_circle
)

old_skater_metrics <- data.frame(
    variable = c(
        "hclust_depth_sd", "hclust_convhull_compact", "hclust_circle_compact"
    ),
    value = c(
        2.5, old_skater_convhull_slope, old_skater_circle_slope
    )
)

# mst_hclust_pixels <- cluster_reef_pixels(
#     pixel_data,
#     habitat_clustering_function = function(x) {
#         mst <- prepare_mst(x)
#         constrained_hclust(x, as_edgelist(mst), alpha = 0.9, beta = -0.85)
#     }
# )
# mst_hclust_pixels$depth <- abs(exact_extract(bathy_raster, mst_hclust_pixels, "mean"))
# mst_hclust_sites <- clustered_pixels_to_polygons(mst_hclust_pixels) %>%
#     mutate(habitat = extract_habitat(site_id)) %>%
#     mutate(depth = abs(exact_extract(bathy_raster, ., "mean")))

# mst_hclust_compact <- st_transform(mst_hclust_sites, crs = 3112) %>%
#     mutate(conv_hull = st_convex_hull(geometry), minbound_circle = st_minimum_bounding_circle(geometry)) %>%
#     mutate(conv_hull_area = st_area(conv_hull), circle_area = st_area(minbound_circle), site_area = st_area(geometry)) %>%
#     mutate(prop_conv_hull = as.numeric(site_area / conv_hull_area), prop_circle = as.numeric(site_area / circle_area))

# mst_hclust_convhull_slope <- mean(
#     mst_hclust_compact[mst_hclust_compact$habitat == "22", ]$prop_conv_hull
# )
# mst_hclust_circle_slope <- mean(
#     mst_hclust_compact[mst_hclust_compact$habitat == "22", ]$prop_circle
# )

params_long <- pivot_longer(params[, c("alpha", "beta", "hclust_depth_sd", "hclust_convhull_compact", "hclust_circle_compact")], c(hclust_depth_sd, hclust_convhull_compact, hclust_circle_compact),names_to="variable", values_to = "value")
ggplot() + geom_line(data = params_long, aes(x = as.numeric(alpha), y = value, color=beta)) +
    geom_hline(data = old_skater_metrics, aes(yintercept = value)) +
    facet_wrap(~variable, scales="free_y")
