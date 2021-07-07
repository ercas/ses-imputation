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
library(terra)
library(sf)

road_buffer_m <- 300
impervious_cutoff <- 5
max_cores <- 30

# State FIPS code to produce crosswalk for. Could theoretically do the entire US
# at once but that would cause long periods of unresponsiveness... also maybe
# can be parallelized this way?
#
# This will be overwritten if a state FIPS code is provided as a CLI arg
#current_statefp <- "25"
current_statefp <- "04"

# Very important! We will be generating very large (~1GB) temporary files -
# don't want to use up all of /tmp/ (default location)
terraOptions(tempdir = "temp/")

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

loadRaster <- function(...) {
  rast(...)
}

# Load data ---------------------------------------------------------------
# We will be transforming everything to the NLCD CRS to ease interoperability
# and ease great-circle distance calculations e.g. buffering

# Detect and validate CLLI FIPS code, if any
args <- commandArgs(trailingOnly=TRUE)
if (length(args) > 0) {
  current_statefp <- args[1]
  if (!file.exists(sprintf("/media/qnap3/ShapeFiles/Polygons/TIGER2000/blocks/tl_2010_%s_tabblock00.shp", current_statefp))) {
    stop(sprintf("no data for FIPS code %s found in /media/qnap3/ShapeFiles/Polygons/TIGER2000/blocks/", current_statefp))
  }
}
message(sprintf("Processing: state with FIPS code %s", current_statefp))

# 2010 Decennial Census data from NHGIS
nhgis_2000_block_pops <- fread('/media/qnap3/Covariates/nhgis/output/2000_block.csv.gz',
                               select = c("GEOID", "population")
                               )[, .(geoid = sprintf("%015.0f", as.numeric(GEOID)),
                                     population = population)]

# 2001 MRLC NLCD land cover
nlcd_land_cover_2001 <- rast("/media/qnap4/MRLC NLCD/landcover/nlcd_2001_land_cover_l48_20210604.img")
nlcd_crs <- crs(nlcd_land_cover_2001)

# 2001 MRLC NLCD impervious surface
nlcd_impervious_2001 <- rast("/media/qnap4/MRLC NLCD/impervious/nlcd_2001_impervious_l48_20210604.img")

# 2000 TIGER geographies to be crosswalked
tiger_2000 <- st_read("/media/qnap3/ShapeFiles/Polygons/TIGER2000/zcta5/tl_2010_us_zcta500.shp") %>%
  st_transform(nlcd_crs)

# 2010 TIGER geographies to use as reference
tiger_2010 <- st_read("/media/qnap3/ShapeFiles/Polygons/TIGER2010/zcta5/tl_2010_us_zcta510.shp") %>%
  st_transform(nlcd_crs)

## State-specific data ----

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
reference_raster <- crop(
  nlcd_land_cover_2001,
  ext(st_buffer(tiger_tracts_2000, road_buffer_m))
)

# We need to create a new variable here called rasterized_id because the raster
# package has some issues saving GEOIDs to GeoTIFFs
blocks_to_rasterize <- load_cached_data(
  sprintf("output/%s_blocks_2000.gpkg", current_statefp),
  function() {
    tiger_blocks_2000 %>%
      transmute(geoid = BLKIDFP00) %>%
      left_join(nhgis_2000_block_pops, on = "geoid") %>%
      mutate(rasterized_id = 1:nrow(tiger_blocks_2000)) %>%
      return()
  },
  save_function = st_write,
  load_function = st_read
)

## Rasterize blocks ----

blocks_rasterized <- load_cached_data(
  sprintf("output/%s_blocks_2000.tif", current_statefp),
  function() {
    message("Rasterizing blocks")
   blocks_to_rasterize_vect <- vect(blocks_to_rasterize)
    return(rasterize(
        blocks_to_rasterize_vect,
        crop(reference_raster, blocks_to_rasterize_vect),
        "rasterized_id"
    ))
  },
  save_function = writeRaster,
  load_function = loadRaster
)

## Buffer and rasterize roads ----
# Here we will be looping over tracts because we can afford to process larger
# areas and move a little faster

roads_rasterized <- load_cached_data(
  sprintf("output/%s_roads_buffered.tif", current_statefp),
  function() {
    message("Buffering roads")
    roads_to_rasterize <- tiger_roads_2010 %>%
        pull(geometry) %>%
        st_buffer(road_buffer_m) %>%
        st_union() %>%
        st_as_sf() %>%
        vect()

    message("Rasterizing roads")
    return(rasterize(
        roads_to_rasterize,
        crop(reference_raster, roads_to_rasterize)
    ))
  },
  save_function = writeRaster,
  load_function = loadRaster
)

## Calculation of inhabited areas ----

inhabited_areas <- load_cached_data(
  sprintf("output/%s_inhabited_areas_2000.tif", current_statefp),
  function() {
    # Start with the road buffer
    message("Adjusting road buffer extents")
    result <- extend(crop(roads_rasterized, blocks_rasterized), blocks_rasterized)
    
    # Remove all areas with <5% impervious surface
    message("Removing <5% impervious surface")
    nlcd_impervious_2001_cropped <- crop(nlcd_impervious_2001, blocks_rasterized)
    stopifnot(dim(result) == dim(nlcd_impervious_2001_cropped))
    result[nlcd_impervious_2001_cropped < 5] <- NA
    
    # Remove all water bodies
    # https://www.mrlc.gov/sites/default/files/NLCD_Colour_Classification_Update.jpg
    message("Removing water bodies")
    nlcd_land_cover_2001_cropped <- crop(nlcd_land_cover_2001, blocks_rasterized)
    stopifnot(dim(result) == dim(nlcd_land_cover_2001_cropped))
    result[nlcd_impervious_2001_cropped == 11] <- NA
    
    # Subset rasterized blocks to this
    message("Masking rasterized blocks")
    result <- mask(blocks_rasterized, result)
    
    return(result)
  },
  save_function = writeRaster,
  load_function = loadRaster
)

## BD interpolation of population counts ----
  
bd_raster <- load_cached_data(
  sprintf("output/%s_population_2000.tif", current_statefp),
  function() {
    # Start with the inhabited areas from earlier - read from disk so we aren't
    # manipulating the original values
    result <- loadRaster(
      sprintf("output/%s_inhabited_areas_2000.tif", current_statefp)
    )
    
    # Divide up block populations evenly across all relevant 30-meter grid cells
    message("Spreading population across grid cells")
    cell_pops <- table(result[]) %>%
      as.data.frame() %>%
      transmute(
        rasterized_id = as.numeric(as.character(Var1)),
        cells = Freq
      ) %>%
      left_join(st_drop_geometry(blocks_to_rasterize), by = "rasterized_id") %>%
      mutate(cell_population = population / cells)
    
    result[] <- left_join(
      data.frame(rasterized_id = result[]),
      cell_pops,
      by = "rasterized_id"
    )$cell_population
    
    return(result)
  },
  save_function = writeRaster,
  load_function = loadRaster
)

# Diagnostics - unreachable code so doesn't run in non-interactive mode
if (FALSE) {
  missing_geoids <- setdiff(blocks_to_rasterize$geoid, cell_pops$geoid)
  length(missing_geoids)
  blocks_to_rasterize %>%
    filter(geoid %in% missing_geoids) %>%
    pull(population) %>%
    #summary()
    sum(na.rm = TRUE)
  
  blocks_to_rasterize %>%
    filter(geoid %in% missing_geoids, population == 0) %>%
    nrow()
}

# Target-Density Weighting (TDW) weights calculation ----------------------

# BD-TDW crosswalk --------------------------------------------------------

# WIP / scratchpad --------------------------------------------------------

