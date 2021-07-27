# Target-Density Weighting (modified) of 2000 -> 2010 Census geographies
#
# General formula: \hat{y}_t = \sum_s{\frac{\left(\frac{A_st}{A_t}\right)z_t}{\sum_\tau\left(\frac{A_{s\tau}}{A_\tau}\right)}y_s} where
# * s = The "source zone", i.e. the data set where the original geographies came
#       from, i.e. the 2000 Census.
# * t = The "target zone", i.e. the data set where the geographies that we want
#       to interpolate to came from, i.e. the 2010 Census.
# * \hat{y}_t = The estimated value of the given Census variable in the target
#       zone, i.e. the value according to 2010 Census geographies.
# * y_s = The value of the given Census variable in the source zone, i.e. the
#       original value in the original geographies (2000 Census).
# * z_t = The value of the given Census variable in the target zone, i.e. the
#       value in a 2010 Census geography.
# * z_\tau = Similar to z_t, but idexing target zone intersections with \tau.
# * A_st = The area of a single intersection between s and t, indexing source
#       zone intersections. There may be multiple, e.g. if a 2000 geography
#       intersects 3 different 2010 geographies.
# * A_s = The area of a geography from the source zone, i.e. a 2000 geography.
# * A_s\tau = Similar to A_st, but indexing target zone intersections instead.
#       There may be multiple, e.g. if a 2010 geography intersects 3 different
#       2000 geographies.
# * A_t = The area of a geography from the target zone, i.e. a 2010 geography.
#
# For reference, the simple areal-weighting formula is \hat{y}_t = \sum_s{\frac{A_{st}}{A_s}y_s}
#
# Primary reference material:
#
# 1. Ruther, M., Leyk, S., & Buttenfield, B. P. (2015). Comparing the effects of
# an NLCD-derived dasymetric refinement on estimation accuracies for multiple
# areal interpolation methods. GIScience & Remote Sensing, 52(2), 158–178.
# https://doi.org/10.1080/15481603.2015.1018856
# 2. Schroeder, J. P. 2007. “Target-Density Weighting Interpolation and
# Uncertainty Evaluation for Temporal Analysis of Census Data.” Geographical
# Analysis 39: 311–335. http://dx.doi.org/10.1111/j.1538-4632.2007.00706.x

#library(doSNOW)
#library(parallel)
library(progress)

library(data.table)
library(dplyr)
library(exactextractr)
library(hash) # For hash tables of heterogeneous type and/or value length
library(hashmap) # For scalar hash tables of homogeneous type
library(terra)
library(sf)

bd_output_directory <- "output/crosswalks/"

terraOptions(tempdir = "temp/")

dir.create(bd_output_directory, showWarnings = FALSE)

load_cached_data <- function(path_rds,
                             f,
                             ...,
                             save_function = saveRDS,
                             load_function = readRDS) {
  if (file.exists(path_rds)) {
    message(sprintf("Loading cached data from %s", path_rds))
    return(load_function(path_rds))
  } else {
    message("No cached data found; running supplied function")
    result <- f(...)
    message(sprintf("Caching data to %s", path_rds))
    save_function(result, path_rds)
    return(result)
  }
}

# Load data ---------------------------------------------------------------

# 2001 MRLC NLCD land cover
nlcd_land_cover_2001 <- rast("/media/qnap4/MRLC NLCD/landcover/nlcd_2001_land_cover_l48_20210604.img")
nlcd_crs <- crs(nlcd_land_cover_2001)

# Decennial Census populations
pops_2000 <- fread("/media/qnap3/Covariates/nhgis/output/2000_zcta.csv.gz"
                   )[, .(ZCTA5CE00 = sprintf("%05.f", as.numeric(GEOID)),
                         pop_sum = population)]
pops_2010 <- fread("/media/qnap3/Covariates/nhgis/output/2010_zcta.csv.gz"
                   )[, .(ZCTA5CE10 = sprintf("%05.f", as.numeric(GEOID)),
                         pop_sum = population)]

# 2000 TIGER geographies to be crosswalked
tiger_2000 <- load_cached_data(
  "output/zctas/zcta5_2000.gpkg",
  function() {
    result <- st_read("/media/qnap3/ShapeFiles/Polygons/TIGER2000/zcta5/tl_2010_us_zcta500.shp") %>%
      st_make_valid()
    result %>%
      st_transform(nlcd_crs) %>%
      st_write("output/zctas/zcta5_2000_nlcd.gpkg")
    return(result)
  },
  save_function = st_write,
  load_function = st_read
)
tiger_2000_nlcd <- load_cached_data(
  "output/zctas/zcta5_2000_nlcd.gpkg",
  function() {
    tiger_2000 %>%
      st_transform(nlcd_crs) %>%
      return()
  },
  save_function = st_write,
  load_function = st_read
)

# 2010 TIGER geographies to use as reference
tiger_2010 <- load_cached_data(
  "output/zctas/zcta5_2010.gpkg",
  function() {
    st_read("/media/qnap3/ShapeFiles/Polygons/TIGER2010/zcta5/tl_2010_us_zcta510.shp") %>%
      st_make_valid() %>%
      return()
  },
  save_function = st_write,
  load_function = st_read
)

tiger_2010_nlcd <- load_cached_data(
  "output/zctas/zcta5_2010_nlcd.gpkg",
  function() {
    tiger_2010 %>%
      st_transform(nlcd_crs) %>%
      return()
  },
  save_function = st_write,
  load_function = st_read
)

# BD-interpolated population from stage 1
bd_pops <- rast("output/bd-pops-conus.tif")
bd_popdens <- rast("output/bd-popdens-conus.tif")

## Calculate intersections ----

tiger_intersections <- load_cached_data(
  "output/zctas/zcta5_2000_2010_intersections.gpkg",
  function() {
    # PostGIS equivalent is much faster - see stage2alt_intersections.sql
    # Attempted parallel implementation with %dopar%; stil slower than PostGIS
    
    st_intersection(tiger_2000["ZCTA5CE00"], tiger_2010["ZCTA5CE10"]) %>%
      st_buffer(0) %>%
      st_cast("MULTIPOLYGON") %>%
      return()
  },
  save_function = st_write,
  load_function = st_read
)
tiger_intersections_nlcd <- load_cached_data(
  "output/zctas/zcta5_2000_2010_intersections_nlcd.gpkg",
  function() {
    tiger_intersections %>%
      st_transform(nlcd_crs) %>%
      return()
  },
  save_function = st_write,
  load_function = st_read
)
tiger_intersections_nlcd$ID <- as.character(tiger_intersections_nlcd$ID)

## Calculate BD-refined areas (inhabited zones) ----

# Bash (from output/ directory):
# $ exactextract --raster pop:bd-pops-conus.tif \
#       --polygons zctas/zcta5_2000_nlcd.gpkg \
#       --fid ZCTA5CE00 \
#       --stat "area_meters=900*count(pop)" \
#       --output zcta-bd-areas/zcta5_2000.csv \
#       --progress | tqdm --total 32038
bd_areas_2000 <- load_cached_data(
  "output/zcta-bd-areas/zcta5_2000.csv",
  function() {
    return(data.table(
      ZCTA5CE00 = tiger_2000_nlcd$ZCTA5CE00,
      area_meters = 900 * exact_extract(bd_pops_raster, tiger_2000_nlcd, "count")
    ))
  },
  save_function = fwrite,
  load_function = fread
)
bd_areas_2000[, ZCTA5CE00 := sprintf("%05.f", as.numeric(ZCTA5CE00))]

# Bash (from output/ directory):
# $ exactextract --raster pop:bd-pops-conus.tif \
#       --polygons zctas/zcta5_2010_nlcd.gpkg \
#       --fid ZCTA5CE10 \
#       --stat "area_meters=900*count(pop)" \
#       --output zcta-bd-areas/zcta5_2010.csv \
#       --progress | tqdm --total 32038
bd_areas_2010 <- load_cached_data(
  "output/zcta-bd-areas/zcta5_2010.csv",
  function() {
    return(data.table(
      ZCTA5CE10 = tiger_2010_nlcd$ZCTA5CE10,
      area_meters = 900 * exact_extract(bd_pops_raster, tiger_2010_nlcd, "count")
    ))
  },
  save_function = fwrite,
  load_function = fread
)
bd_areas_2010[, ZCTA5CE10 := sprintf("%05.f", as.numeric(ZCTA5CE10))]

# Bash (from output/ directory):
# $ exactextract --raster pop:bd-pops-conus.tif \
#       --polygons zctas/zcta5_2000_2010_intersections_nlcd.gpkg \
#       --fid ID \
#       --stat "area_meters=900*count(pop)" \
#       --output zcta-bd-areas/zcta5_2000_2010_intersections.csv \
#       --progress | tqdm --total 133205
bd_areas_2000_2010 <- load_cached_data(
  "output/zcta-bd-areas/zcta5_2000_2010_intersections.csv",
  function() {
    return(data.table(
      ID = tiger_intersections_nlcd$ID,
      area_meters = 900 * exact_extract(bd_pops_raster, tiger_intersections_nlcd, "count")
    ))
  },
  save_function = fwrite,
  load_function = fread
)
bd_areas_2000_2010[, ID := as.character(ID)]

## Calculate BD-refined areas, lenient definition ----

bd_areas_2000_lenient <- load_cached_data(
  "output/zcta-bd-areas/zcta5_2000.csv",
  function() {
    return(data.table(
      ZCTA5CE00 = tiger_2000_nlcd$ZCTA5CE00,
      area_meters = 900 * exact_extract(bd_pops_raster, tiger_2000_nlcd, "count")
    ))
  },
  save_function = fwrite,
  load_function = fread
)
bd_areas_2000_lenient[, ZCTA5CE00 := sprintf("%05.f", as.numeric(ZCTA5CE00))]

bd_areas_2010_lenient <- load_cached_data(
  "output/zcta-bd-areas/zcta5_2010.csv",
  function() {
    return(data.table(
      ZCTA5CE10 = tiger_2010_nlcd$ZCTA5CE10,
      area_meters = 900 * exact_extract(bd_pops_raster, tiger_2010_nlcd, "count")
    ))
  },
  save_function = fwrite,
  load_function = fread
)
bd_areas_2010_lenient[, ZCTA5CE10 := sprintf("%05.f", as.numeric(ZCTA5CE10))]

bd_areas_2000_2010 <- load_cached_data(
  "output/zcta-bd-areas/zcta5_2000_2010_intersections.csv",
  function() {
    return(data.table(
      ID = tiger_intersections_nlcd$ID,
      area_meters = 900 * exact_extract(bd_pops_raster, tiger_intersections_nlcd, "count")
    ))
  },
  save_function = fwrite,
  load_function = fread
)
bd_areas_2000_2010_lenient[, ID := as.character(ID)]

## Calculate BD-interpolated populations ----
# Deprecated section - intersection population is not needed, so we can just
# pull counts directly from the decennial Census since all population data will
# be from the same universe

# Bash (from output/ directory):
# $ exactextract --raster pop:bd-pops-conus.tif \
#       --polygons zctas/zcta5_2000_nlcd.gpkg \
#       --fid ZCTA5CE00 \
#       --stat "sum(pop)" \
#       --output zcta-bd-pops/zcta5_2000.csv \
#       --progress | tqdm --total 32038
# pops_2000 <- load_cached_data(
#   "output/zcta-bd-pops/zcta5_2000.csv",
#   function() {
#     return(data.table(
#       ZCTA5CE00 = tiger_2000_nlcd$ZCTA5CE00,
#       pop_sum = exact_extract(bd_pops_raster, tiger_2000_nlcd, "sum")
#     ))
#   },
#   save_function = fwrite,
#   load_function = fread
# )
# pops_2000[, ZCTA5CE00 := sprintf("%05.f", as.numeric(ZCTA5CE00))]

# Bash (from output/ directory):
# $ exactextract --raster pop:bd-pops-conus.tif \
#       --polygons zctas/zcta5_2010_nlcd.gpkg \
#       --fid ZCTA5CE10 \
#       --stat "sum(pop)" \
#       --output zcta-bd-pops/zcta5_2010.csv \
#       --progress | tqdm --total 33120
# pops_2010 <- load_cached_data(
#   "output/zcta-bd-pops/zcta5_2010.csv",
#   function() {
#     return(data.table(
#       ZCTA5CE10 = tiger_2010_nlcd$ZCTA5CE10,
#       pop_sum = exact_extract(bd_pops_raster, tiger_2010_nlcd, "sum")
#     ))
#   },
#   save_function = fwrite,
#   load_function = fread
# )
# pops_2010[, ZCTA5CE10 := sprintf("%05.f", as.numeric(ZCTA5CE10))]

# Bash (from output/ directory):
# $ exactextract --raster pop:bd-pops-conus.tif \
#       --polygons zctas/zcta5_2000_2010_intersections_nlcd.gpkg \
#       --fid ID \
#       --stat "sum(pop)" \
#       --output zcta-bd-pops/zcta5_2000_2010_intersections.csv \
#       --progress | tqdm --total 133205
# pops_2000_2010 <- load_cached_data(
#   "output/zcta-bd-pops/zcta5_2000_2010_intersections.csv",
#   function() {
#     return(data.table(
#       ID = as.character(tiger_intersections_nlcd$ID),
#       pop_sum = exact_extract(bd_pops_raster, tiger_intersections_nlcd, "sum")
#     ))
#   },
#   save_function = fwrite,
#   load_function = fread
# )
# pops_2000_2010 <- pops_2000_2010[st_drop_geometry(tiger_intersections_nlcd), on = "ID"
#                                  ][, `:=` (ID = as.character(ID),
#                                            ZCTA5CE00 = sprintf("%05.f", as.numeric(ZCTA5CE00)),
#                                            ZCTA5CE10 = sprintf("%05.f", as.numeric(ZCTA5CE10)))]

# Target-Density Weighting ------------------------------------------------
# There's definitely a more efficient way to do this

intersections <- as.data.table(st_drop_geometry(tiger_intersections)
                               )[, `:=` (ID = as.character(ID),
                                         ZCTA5CE00 = as.character(ZCTA5CE00),
                                         ZCTA5CE10 = as.character(ZCTA5CE10))
                                 ][bd_areas_2000_2010, on = "ID"]

## Create hash tables to speed up queries ----

# Map of 2000 ZCTAs -> 2000 ZCTA areas
bd_areas_2000_hash <- hashmap(bd_areas_2000$ZCTA5CE00, bd_areas_2000$area_meters)

# Map of 2010 ZCTAs -> 2010 ZCTA areas
bd_areas_2010_hash <- hashmap(bd_areas_2010$ZCTA5CE10, bd_areas_2010$area_meters)

# Map of intersection ID -> 2000-2010 ZCTA intersection areas
bd_areas_2000_2010_hash <- hashmap(bd_areas_2000_2010$ID, bd_areas_2000_2010$area_meters)

# Map of 2000 ZCTAs -> 2000 BD-interpolated populations (extraneous?)
pops_2000_hash <- hashmap(pops_2000$ZCTA5CE00, pops_2000$pop_sum)

# Map of 2010 ZCTAs -> 2010 BD-interpolated populations (extraneous?)
pops_2010_hash <- hashmap(pops_2010$ZCTA5CE10, pops_2010$pop_sum)

# Map of 2010 ZCTAs -> vector of intersecting 2000 ZCTAs
zctas_2010_to_2000_hash <- hash()
for (zcta in tiger_2010$ZCTA5CE10) {
  zctas_2010_to_2000_hash[[zcta]] <- as.character(
    intersections[ZCTA5CE10 == zcta]$ZCTA5CE00
  )
}

# Map of 2000 ZCTAs -> vector of intersecting 2010 ZCTAs
zctas_2000_to_2010_hash <- hash()
for (zcta in tiger_2000$ZCTA5CE00) {
  zctas_2000_to_2010_hash[[zcta]] <- as.character(
    intersections[ZCTA5CE00 == zcta]$ZCTA5CE10
  )
}

# Map of 2000 ZCTAs -> vector of intersecting 2010 ZCTAs
intersection_id_2000_to_2010_hash <- hash()
for (zcta in tiger_2000$ZCTA5CE00) {
  nested_hash <- hash()
  temp <- intersections[ZCTA5CE00 == zcta]
  for (i in 1:nrow(temp)) {
    row <- temp[i,]
    nested_hash[[row$ZCTA5CE10]] <- row$ID
  }
  intersection_id_2000_to_2010_hash[[zcta]] <- nested_hash
}

# Make sure hash tables are working
stopifnot(sum(is.na(lapply(bd_areas_2000$ZCTA5CE00, function(x) bd_areas_2000_hash[[x]]))) == sum(is.na(bd_areas_2000$area_meters)))
stopifnot(sum(is.na(lapply(bd_areas_2010$ZCTA5CE10, function(x) bd_areas_2010_hash[[x]]))) == sum(is.na(bd_areas_2010$area_meters)))
stopifnot(sum(is.na(lapply(bd_areas_2000_2010$ID, function(x) bd_areas_2000_2010_hash[[x]]))) == sum(is.na(bd_areas_2000_2010$area_meters)))
stopifnot(sum(is.na(lapply(pops_2000$ZCTA5CE00, function(x) pops_2000_hash[[x]]))) == sum(is.na(pops_2000$pop_sum)))
stopifnot(sum(is.na(lapply(pops_2010$ZCTA5CE10, function(x) pops_2010_hash[[x]]))) == sum(is.na(pops_2010$pop_sum)))
stopifnot(sum(is.na(lapply(intersections$ZCTA5CE00, function(x) zctas_2000_to_2010_hash[[x]]))) == 0)
stopifnot(sum(is.na(lapply(intersections$ZCTA5CE10, function(x) zctas_2010_to_2000_hash[[x]]))) == 0)

## Calculation of weights ----

calculate_tdw <- function(areas_2000_hash = bd_areas_2000_hash,
                          areas_2010_hash = bd_areas_2010_hash,
                          areas_2000_2010_hash = bd_areas_2000_2010_hash) {
  bar <- progress_bar$new(
    "Calculating target-density weights :current/:total (:percent) [:bar] eta :eta",
    total = nrow(tiger_2010)
  )
  result <- lapply(
    as.character(tiger_2010$ZCTA5CE10),
    function(target_zcta) {
      # A_t
      target_area <- areas_2010_hash[[target_zcta]]
      
      # Z_t
      target_pop <- pops_2010_hash[[target_zcta]]
      
      # s (index)
      source_zctas <- zctas_2010_to_2000_hash[[target_zcta]]
      
      # sum_s(A_s,t / A_t * Z_t / sum_tau(A_s,tau / A_tau * Z_tau))
      weights <- sapply( # Sum over intersecting source ZCTAs (index s)
        source_zctas,
        function(source_zcta) {
          # A_s,t
          #intersection_area <- intersections[ZCTA5CE10 == target_zcta & ZCTA5CE00 == source_zcta]$area_meters
          intersection_id <- intersection_id_2000_to_2010_hash[[source_zcta]][[target_zcta]]
          intersection_area <- areas_2000_2010_hash[[intersection_id]]
          
          # A_s,t / A_t * Z_t
          numerator <- intersection_area / target_area * target_pop
          
          denominator <- sum(sapply(
            # tau (index)
            zctas_2000_to_2010_hash[[source_zcta]],
            function(target_zcta_2) { # Sum over intersecting target ZCTAs (index tau)
              # A_s,tau
              #intersection_area_2 <- intersections[ZCTA5CE10 == target_zcta_2 & ZCTA5CE00 == source_zcta]$area_meters
              intersection_id_2 <- intersection_id_2000_to_2010_hash[[source_zcta]][[target_zcta_2]]
              intersection_area_2 <- areas_2000_2010_hash[[intersection_id_2]]
              
              # A_tau
              target_area_2 <- areas_2010_hash[[target_zcta_2]]
              
              # Z_tau
              target_pop_2 <- pops_2010_hash[[target_zcta_2]]
              
              # A_s,tau / A_tau * Z_tau
              return(intersection_area_2 / target_area_2 * target_pop_2)
            }
          ))
          
          # A_s,t / A_t * Z_t / sum_tau(A_s,tau / A_tau * Z_tau)
          return(numerator / denominator)
        }
      )
      
      bar$tick()
      
      result <- data.table(
        ZCTA5CE10 = target_zcta,
        ZCTA5CE00 = source_zctas,
        tdw = weights
      )
      if (all(is.na(result$ZCTA5CE00))) {
        return(NA)
      } else {
        return(result)
      }
    }
  )
  result <- bind_rows(result[!is.na(result)])
  return(result)
}

# Pass 1: fully bd-refined
tdw <- load_cached_data(
  "output/tdw.csv",
  function() {
    calculate_tdw(
      areas_2000_hash = bd_areas_2000_hash,
      areas_2010_hash = bd_areas_2010_hash,
      areas_2000_2010_hash = bd_areas_2000_2010_hash,
    )
  }
  save_function = fwrite,
  load_function = fread
)
tdw[, `:=` (ZCTA5CE10 = sprintf("%05.f", as.numeric(ZCTA5CE10)),
            ZCTA5CE00 = sprintf("%05.f", as.numeric(ZCTA5CE00)))]

# Pass 2:

# Diagnostics -------------------------------------------------------------

library(ggplot2)
library(rasterVis) # For gplot
library(viridis)

tiger_2000_dt <- as.data.table(tiger_2000)
tiger_2000_nlcd_dt <- as.data.table(tiger_2000_nlcd)

plot_tdw <- function(debug_zcta) {
  # sf of the target ZCTA
  debug_target <- filter(tiger_2010, ZCTA5CE10 == debug_zcta)
  print(debug_target)
  
  # sf of the source ZCTAs with their TDWs
  debug_tdw_sf <- st_as_sf(
    tiger_2000_dt[tdw[ZCTA5CE10 == debug_zcta], on = "ZCTA5CE00"]
  )
  
  ggplot() +
    geom_sf(
      aes(fill = log(tdw)),
      data = debug_tdw_sf,
      lwd = NA
    ) +
    geom_sf(
      data = debug_target,
      fill = NA, color = "red", lwd = 1.5,
    ) +
    geom_sf_text(
      aes(label = ZCTA5CE00),
      data = debug_tdw_sf,
      color = "white"
    ) +
    labs(
      title = sprintf("Target-density weights for %s", debug_zcta),
      subtitle = sprintf("2000 ZCTA GEOIDs (%s) shown in white text", nrow(debug_tdw_sf)),
      x = "Longitude",
      y = "Latitude"
    ) %>%
    return()
}

plot_bd_popdens <- function(debug_zcta) {
  debug_target_nlcd <- filter(tiger_2010_nlcd, ZCTA5CE10 == debug_zcta)
  
  debug_tdw_sf_nlcd <- st_as_sf(
    tiger_2000_nlcd_dt[tdw[ZCTA5CE10 == debug_zcta], on = "ZCTA5CE00"]
  )
  
  # For some reason, st_bbox doesn't give us the correct extents - manual calc
  bounds <- as.data.table(do.call(rbind, lapply(debug_tdw_sf_nlcd$geom, st_bbox)))
  extents <- ext(min(bounds$xmin), max(bounds$xmax), min(bounds$ymin), max(bounds$ymax))
  
  gplot(crop(bd_popdens, extents)) +
    geom_raster(aes(fill = value)) +
    scale_fill_viridis(na.value = NA) +
    geom_sf(
      data = debug_tdw_sf_nlcd,
      fill = NA, color = "orange", lwd = 0.5,
      inherit.aes = FALSE
    ) +
    geom_sf(
      data = debug_target_nlcd,
      fill = NA, color = "red", lwd = 1,
      inherit.aes = FALSE
    ) +
    labs(
      title = sprintf("BD-interp. pop. dens. near %s", debug_zcta),
      subtitle = sprintf("2000 ZCTAs (%s) shown in orange", nrow(debug_tdw_sf_nlcd)),
      x = "Longitude",
      y = "Latitude"
    ) %>%
    return()
}

# Unreachable code so doesn't run in non-interactive mode
if (FALSE) {
  # Harvard T.H. Chan
  plot_tdw("02115")
  plot_bd_popdens("02115")
  
  # Random
  plot_tdw("76034")
  plot_bd_popdens("76034")
  
  # Number of 2000 ZCTAs per 2010 ZCTA
  tdw[, .N, by = "ZCTA5CE10"][order(N)]
  plot_tdw("15767")
  plot_bd_popdens("15767")
  
  # Most even spread of weights
  tdw %>%
    group_by(ZCTA5CE10) %>%
    summarize(
      n = n(),
      median_tdw = median(tdw, na.rm = TRUE),
      mean_tdw = mean(tdw, na.rm = TRUE)
    ) %>%
    filter(
      n > 1,
      !is.na(median_tdw)
    ) %>%
    arrange(-median_tdw) %>%
    View()
  plot_tdw("97434")
  plot_bd_popdens("97434")
  tdw[ZCTA5CE10 == "97434"]
  
  plot_tdw("25661")
  plot_bd_popdens("25661")
  sum(pops_2000[tdw[ZCTA5CE10 == "25661"], on = "ZCTA5CE00"]$pop_sum)
  sum(pops_2000[tdw[ZCTA5CE10 == "25661"], on = "ZCTA5CE00"][is.na(tdw)]$pop_sum)
  sum(pops_2000[tdw[is.na(tdw)], on = "ZCTA5CE00"]$pop_sum) - pop_hi_pr
  sum(pops_2000$pop_sum) - pop_hi_pr
  pop_hi_pr <- sum(pops_2000[bd_areas_2000[is.na(bd_area_meters)], on = "ZCTA5CE00"]$pop_sum)
  
  plot_tdw("")
  plot_bd_popdens("")
  
  plot_tdw(sample(tiger_2010$ZCTA5CE10, 1))
}

