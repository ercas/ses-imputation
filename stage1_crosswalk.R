# General approach: Binary dasymmetric interpolation
#
# This is a simplification of the NHGIS approach, which uses both 
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
#
# Primary reference material:
# 1. Ruther, M., Leyk, S., & Buttenfield, B. P. (2015). "Comparing the effects
# of an NLCD-derived dasymetric refinement on estimation accuracies for multiple
# areal interpolation methods." GIScience & Remote Sensing 52(2), 158-178.
# http://dx.doi.org/10.1080/15481603.2015.1018856
# 2. Schroeder, J. P. 2007. “Target-Density Weighting Interpolation and
# Uncertainty Evaluation for Temporal Analysis of Census Data.” Geographical
# Analysis 39: 311–335. http://dx.doi.org/10.1111/j.1538-4632.2007.00706.x
#
# Contact: Edgar Castro <edgar_castro@g.harvard.edu>

library(doSNOW)
library(parallel)
library(progress)

library(data.table)
library(dplyr)
library(fasterize)
library(ggplot2)
library(raster)
library(sf)

road_buffer_m <- 300
impervious_cutoff <- 5
max_cores <- 30

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
# We will be transforming everything to the NLCD CRS to ease interoperability
# and ease great-circle distance calculations e.g. buffering

# 2010 Decennial Census data from NHGIS
nhgis_2000_block_pops <- fread('/media/qnap3/Covariates/nhgis/output/2000_block.csv.gz',
                               select = c("GEOID", "population")
                               )[, .(geoid = sprintf("%015.0f", as.numeric(GEOID)),
                                     population = population)]

# 2001 MRLC NLCD land cover
nlcd_land_cover_2001 <- brick("/media/qnap4/MRLC NLCD/landcover/nlcd_2001_land_cover_l48_20210604.img")
nlcd_crs <- st_crs(nlcd_land_cover_2001)

# 2001 MRLC NLCD impervious surface
nlcd_impervious_2001 <- brick("/media/qnap4/MRLC NLCD/impervious/nlcd_2001_impervious_l48_20210604.img")

# 2000 TIGER geographies to be crosswalked
tiger_2000 <- st_read("/media/qnap3/ShapeFiles/Polygons/TIGER2000/zcta5/tl_2010_us_zcta500.shp") %>%
  st_transform(nlcd_crs)

# 2010 TIGER geographies to use as reference
tiger_2010 <- st_read("/media/qnap3/ShapeFiles/Polygons/TIGER2010/zcta5/tl_2010_us_zcta510.shp") %>%
  st_transform(nlcd_crs)

## State-specific data ----

# State FIPS code to produce crosswalk for. Could theoretically do the entire US
# at once but that would cause long periods of unresponsiveness... also maybe
# can be parallelized this way?
current_statefp <- "25"

# 2000 TIGER tracts
tiger_tracts_2000 <- sprintf("/media/qnap3/ShapeFiles/Polygons/TIGER2000/tracts/tl_2010_%s_tract00.shp", current_statefp) %>%
  st_read() %>%
  st_transform(nlcd_crs)

# 2000 TIGER blocks
tiger_blocks_2000 <- sprintf("/media/qnap3/ShapeFiles/Polygons/TIGER2000/blocks/tl_2010_%s_tabblock00.shp", current_statefp) %>%
  st_read() %>%
  st_transform(nlcd_crs)

# 2010 TIGER blocks
tiger_blocks_2010 <- sprintf("/media/qnap3/ShapeFiles/Polygons/TIGER2010/blocks/tl_2010_%s_tabblock10.shp", current_statefp) %>%
  st_read() %>%
  st_transform(nlcd_crs)

# 2010 TIGER roads
tiger_roads_2010 <- sprintf("/media/qnap3/ShapeFiles/Lines/TIGER2010/roads/tl_2010_%s*_roads.shp", current_statefp) %>%
  Sys.glob() %>%
  lapply(st_read) %>%
  bind_rows() %>%
  st_transform(nlcd_crs)

# Binary dasymetric (BD) interpolation ------------------------------------

# Load the entire reference area into memory for faster processing
reference_raster <- readAll(crop(
  nlcd_land_cover_2001,
  extent(st_buffer(tiger_tracts_2000, road_buffer_m))
))$Layer_1

## Rasterize blocks ----

blocks_to_rasterize <- tiger_blocks_2000 %>%
  transmute(geoid = BLKIDFP00) %>%
  left_join(nhgis_2000_block_pops, on = "geoid") %>%
  mutate(geoid = as.numeric(geoid)) # for raster storage

bar <- progress_bar$new(
  "* Rasterizing geographies :current/:total (:percent) [:bar] eta :eta",
  total = nrow(blocks_to_rasterize)
)
blocks_rasterized <- lapply(
  1:nrow(blocks_to_rasterize),
  function(i) {
    bar$tick()
    block <- blocks_to_rasterize[i,]
    return(fasterize(block, crop(reference_raster, block), "geoid"))
  }
)

message("Merging and cropping rasterized geographies")
blocks_rasterized <- crop(do.call(merge, blocks_rasterized), blocks_to_rasterize)

writeRaster(blocks_rasterized, sprintf("output/%s_blocks_2000.tif", current_statefp))

## Buffer and rasterize roads ----
# Here we will be looping over tracts because we can afford to process larger
# areas and move a little faster

bar <- progress_bar$new(
  "* Rasterizing roads :current/:total (:percent) [:bar] eta :eta",
  total = nrow(tiger_tracts_2000)
)
roads_rasterized <- lapply(
  1:nrow(tiger_tracts_2000),
  function(i) {
    bar$tick()
    tract <- tiger_tracts_2000[i,]
    # tryCatch block because some areas may have no roads at all; leads to
    # Error in fasterize(...): sf geometry must be POLYGON or MULTIPOLYGON
    tryCatch(
      {
        road_buffer <- tiger_roads_2010 %>%
          st_intersection(tiger_tracts_2000[i,]) %>%
          pull(geometry) %>%
          st_buffer(road_buffer_m) %>%
          st_union() %>%
          st_as_sf()
        return(fasterize(road_buffer, crop(reference_raster, road_buffer)))
      },
      error = function(...) return(NA)
    )
  }
)

# May be a good time to check on those null indices before they get omitted and
# confirm that it was due to no roads
roads_rasterized <- na.omit(roads_rasterized)

message("Merging rasterized road buffers")
roads_rasterized <- do.call(merge, roads_rasterized)

writeRaster(roads_rasterized, sprintf("output/%s_roads_buffered.tif", current_statefp))

## BD interpolation of population counts ----

# Start with the road buffer
bd_raster <- extend(crop(roads_rasterized, blocks_rasterized), blocks_rasterized)
plot(bd_raster)

# Remove all areas with <5% impervious surface
nlcd_impervious_2001_cropped <- crop(nlcd_impervious_2001, blocks_rasterized)
stopifnot(compareRaster(bd_raster, nlcd_impervious_2001_cropped))
bd_raster[nlcd_impervious_2001_cropped < 5] <- NA
plot(bd_raster)

# Remove all water bodies
# https://www.mrlc.gov/sites/default/files/NLCD_Colour_Classification_Update.jpg
nlcd_land_cover_2001_cropped <- crop(nlcd_land_cover_2001, blocks_rasterized)
stopifnot(compareRaster(bd_raster, nlcd_land_cover_2001_cropped))
bd_raster[nlcd_impervious_2001_cropped == 11] <- NA
plot(bd_raster)

# Subset rasterized blocks to this
bd_raster <- mask(blocks_rasterized, bd_raster)
plot(bd_raster)

writeRaster(roads_rasterized, sprintf("output/%s_inhabited_areas_2000.tif", current_statefp))

# Divide up block populations evenly across all relevant 30-meter grid cells
cell_pops <- table(bd_raster[]) %>%
  as.data.frame() %>%
  transmute(
    geoid = as.numeric(as.character(Var1)),
    cells = Freq
  ) %>%
  left_join(st_drop_geometry(blocks_to_rasterize)) %>%
  mutate(cell_population = population / cells)
bd_raster[] <- left_join(
  data.frame(geoid = bd_raster[]),
  cell_pops
)$cell_population
plot(bd_raster)

writeRaster(roads_rasterized, sprintf("output/%s_population_2000.tif", current_statefp))

# Target-Density Weighting (TDW) weights calculation ----------------------

# BD-TDW crosswalk --------------------------------------------------------

# WIP / scratchpad --------------------------------------------------------

x1 <- tiger_blocks_2000[1,]
x2 <- fasterize(x1, crop(nlcd_land_cover_2001$Layer_1, x1), "BLKIDFP00")
y1 <- tiger_blocks_2000[2,]
y2 <- fasterize(y1, crop(nlcd_land_cover_2001$Layer_1, y1), "BLKIDFP00")
z1 <- head(tiger_blocks_2000, 100)
z1$geoid <- as.numeric(as.character(z1$BLKIDFP00))
z2 <- fasterize(z1, crop(nlcd_land_cover_2001$Lazer_1, z1), "geoid")
plot(merge(x2, y2))

test_zcta <- tiger_2000[tiger_2000$ZCTA5CE00 == "01960",]

plot(crop(nlcd_impervious_2001, test_zcta))
plot(crop(nlcd_impervious_2001, test_zcta) >= impervious_cutoff)

# Buffer the road geometries
tiger_roads_2010 %>%
  st_intersection(test_zcta) %>%
  pull(geometry) %>%
  st_buffer(road_buffer_m) %>%
  st_union() %>%
  plot()

nlcd_impervious_2001 %>%
  crop(test_zcta) %>%
  mask(test_zcta) %>%
  rasterToPolygons(function(impervious) impervious >= impervious_cutoff) %>%
  st_as_sf() %>%
  st_union() %>%
  plot()

tiger_roads_2010 %>%
  st_intersection(test_zcta) %>%
  pull(geometry) %>%
  st_buffer(road_buffer_m) %>%
  st_union() %>%
  st_difference(
    nlcd_impervious_2001 %>%
      crop(test_zcta) %>%
      rasterToPolygons(function(impervious) impervious >= impervious_cutoff) %>%
      st_as_sf() %>%
      st_intersection(test_zcta) %>%
      st_union()
  ) %>%
  plot()