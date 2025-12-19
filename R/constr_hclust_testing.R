library(ReefPartitionAllenAtlas)

library(sf)
library(terra)
library(h3)
library(tidyverse)
library(igraph)
library(exactextractr)
library(furrr)
library(RColorBrewer)

plan(multisession, workers = 4)

# Set seed for reproducibility
set.seed(123)

ID <- "23055101104" # One Tree Island Reef

reef_polygons <- st_read("../canonical-reefs/output/rrap_canonical_2025-01-10-T13-28-06.gpkg")

bathy_raster <- rast("../../reef_partitioning/Data/required_data/bathy/Mackay-Capricorn/MC_wg84.tif")
habitat_raster <- rast("../../reef_partitioning/Data/required_data/GBR10 GBRMP Geomorphic.tif")

# Must ensure all input data are in the same CRS
reef_polygons <- st_make_valid(st_transform(reef_polygons, crs = sf::st_crs(habitat_raster)))

habitat_categories <- c(15, 22, 14, 21) # reef crest, reef slope, outer reef flat, sheltered reef slope

# Define hexagon resolution and unit for outputting hexagon area (same as function defaults)
hex_resolution <- 12
unit <- "km2"

site_size <- 250 * 250
hex_size <- data.frame(
        Res = c(7:15),
        Size = c(5161293, 737327, 105332, 15047, 2149, 307.09, 43.87, 6.267, 0.895)
    )
pix_per_site <- round(site_size / hex_size$Size[hex_size$Res == hex_resolution])


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
        constrained_hclust(x, as_edgelist(mst), alpha = 0.595)
    }
)
mst_hclust_sites <- clustered_pixels_to_polygons(mst_hclust_pixels)
mst_hclust_sites$depth <- abs(exact_extract(bathy_raster, mst_hclust_sites, "mean"))
ggplot() + geom_sf(data = mst_hclust_sites, aes(fill = depth)) + scale_fill_fermenter(palette = "Greens", direction = 1, breaks = c(2.5, 5, 7.5, 10, 12.5, 15, 17.5, 20))
ggplot() + geom_sf(data = mst_hclust_pixels, aes(color = depth)) + scale_color_fermenter(palette = "Greens", direction = 1, breaks = c(2.5, 5, 7.5, 10, 12.5, 15, 17.5, 20))

mst_hclust_sites$area <- as.numeric(st_area(mst_hclust_sites)) / 1e6

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
        constrained_hclust(x, tri_edges, alpha = 0.595)
    }
)
high_weight_hclust_sites <- clustered_pixels_to_polygons(high_weight_hclust_pixels)

# # Test approach with skater (as in previous functions)
# skater_pixels <- cluster_reef_pixels(
#     pixel_data,
#     habitat_clustering_function = function(x) reef_skater(x)
# )
# saveRDS(skater_pixels, "OTIReef_skater.rds")
# skater_pixels <- readRDS("OTIReef_skater.rds")
# skater_sites <- clustered_pixels_to_polygons(skater_pixels)

extract_habitat <- function(id) {
    one <- str_extract(as.character(id), "_\\d{2}_")
    two <- str_extract(one, "\\d{2}")
}

# skater_compactness <- skater_sites %>%
#     st_transform(crs = 3112) %>%
#     mutate(habitat = extract_habitat(site_id)) %>%
#     mutate(conv_hull = st_convex_hull(geometry), minbound_circle = st_minimum_bounding_circle(geometry)) %>%
#     mutate(conv_hull_area = st_area(conv_hull), circle_area = st_area(minbound_circle), site_area = st_area(geometry)) %>%
#     mutate(prop_conv_hull = site_area / conv_hull_area, prop_circle = site_area / circle_area)

# skater_pixels$hexagons <- h3_to_geo_boundary_sf(skater_pixels$h3_index)
# skater_pixels$depth <- abs(exact_extract(bathy_raster, skater_pixels$hexagons, "mean"))

# skater_variance <- skater_pixels %>%
#     group_by(habitat, site_id) %>%
#     summarise(depth_sd = sd(depth), .groups = "keep") %>%
#     group_by(habitat) %>%
#     summarise(mean_habitat_sd = mean(depth_sd), variance_habitat_sd = var(depth_sd))

# mean_skater_slope_sd <- as.numeric(mean(skater_variance[skater_variance$habitat == 22, ]$mean_habitat_sd))
# mean_skater_convhull_compact <- as.numeric(mean(skater_compactness[skater_compactness$habitat == 22, ]$prop_conv_hull, na.rm=TRUE))
# mean_skater_circle_compact <- as.numeric(mean(skater_compactness[skater_compactness$habitat == 22, ]$prop_circle, na.rm = TRUE))

compare_skater_slopes <- function(param) {
    alpha <- param$alpha
    beta <- param$beta
    method <- param$method

    hclust_pixels <- cluster_reef_pixels(
        pixel_data,
        habitat_clustering_function = function(x) {
            mst <- prepare_mst(x)
            constrained_hclust(
                x, 
                as_edgelist(mst), 
                alpha = alpha, 
                beta = beta, 
                method = method,
                n_clust = round(nrow(x) / pix_per_site)
            )
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

alphas <- seq(0, 1, 0.005)
betas <- -1
methods <- "ward.D2"

params <- expand.grid(alpha = alphas, beta = betas, method = methods)
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
    x[, "method"] <- params_list[[x_ind]]$method
    x
})
param_clusters <- do.call(rbind, param_clusters)

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
        2.51, old_skater_convhull_slope, old_skater_circle_slope
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

params_long <- pivot_longer(params[, c("alpha", "method", "hclust_depth_sd", "hclust_convhull_compact", "hclust_circle_compact")], c(hclust_depth_sd, hclust_convhull_compact, hclust_circle_compact),names_to="variable", values_to = "value")
ggplot() + geom_line(data = params_long, aes(x = as.numeric(alpha), y = value, color=method)) +
    geom_hline(data = old_skater_metrics, aes(yintercept = value)) +
    facet_wrap(~variable, scales="free_y")

summary_vals <- group_by(params, alpha, method) %>%
    mutate(
        sd_diff = abs(2.51 - hclust_depth_sd), 
        convhull_diff = abs(old_skater_convhull_slope - hclust_convhull_compact), 
        circle_diff = abs(old_skater_circle_slope - hclust_circle_compact)
    ) %>%
    mutate(sum_diff = sd_diff + convhull_diff + circle_diff)



# Benchmarking constr.hclust performance across reef sizes once parameters are chosen
bathy_box <- st_as_sfc(st_bbox(bathy_raster))
within <- sapply(reef_polygons$geometry, function(x) st_within(st_make_valid(x), bathy_box))
within <- which(within == 1)

within_ids <- reef_polygons[within, ]$UNIQUE_ID
within_ids <- sample(within_ids, 700)

reef_cluster_time <- function(pix) {
    if (is.null(pix)) {return(NULL)}
    print(unique(pix$UNIQUE_ID))
    habitat_more_20 <- names(table(pix$habitat))[table(pix$habitat) > 20]

    pix <- pix[pix$habitat %in% habitat_more_20, ]

    mst_hclust_pixels <- NULL
    tryCatch(mst_hclust_pixels <- cluster_reef_pixels(
        pix,
        habitat_clustering_function = function(x) {
            mst <- prepare_mst(x)

            start_time <- Sys.time()
            clusters <- constrained_hclust(x, as_edgelist(mst), alpha=0.595, n_clust = max(1, round(nrow(x) / pix_per_site)))
            end_time <- Sys.time()
            clusters$cluster_time <- end_time - start_time

            clusters
        }
    ), error = function(e) {return(NULL)}
    )
    if (is.null(mst_hclust_pixels)) {return(NULL)}
    gc()
    mst_hclust_pixels
}

target_reefs <- reef_polygons[reef_polygons$UNIQUE_ID %in% within_ids, ]
pixel_data <- sapply(split(target_reefs, target_reefs$UNIQUE_ID), function(x) {
    print(unique(x$UNIQUE_ID))
    pixels <- NULL
    tryCatch(pixels <- extract_pixel_points(
            x,
            habitat_raster,
            bathy_raster,
            habitat_categories,
            hex_resolution = hex_resolution,
            unit = unit,
            additional_variable_name = "depth",
            output_epsg = 3112,
            resample_method = "bilinear"
        ), error = function(e) {return(NULL)}
    )
    if (is.null(pixels)) {return(NULL)}
    if (nrow(pixels) < 4) {return(NULL)}

    pixels$UNIQUE_ID <- unique(x$UNIQUE_ID)

    pixels
})

pixel_data <- pixel_data[!sapply(pixel_data, is.null)]
hclust_times <- lapply(
    pixel_data, 
    reef_cluster_time
)

hclust_times <- hclust_times[!sapply(hclust_times, is.null)]

clustering_times <- lapply(hclust_times, function(x) {
    if(nrow(x) == 0) {return()}
    print(unique(x$UNIQUE_ID))
    times <- unique(x[, c("cluster_time", "npixels"), drop = TRUE])

    times
})
clustering_times <- do.call(rbind, clustering_times)

problem_id <- "19130100104"

