# Binary dasymetric interpolation of Census population
#
# This script generates a population surface at the resolution of the MRLC NLCD
# 2001 (30-meter grid) with population interpolated dasymetrically across each
# block group. In short, the process is as follows:
#
# 1. Retrieve data:
#     * Block-level population counts from the 2000 decennial Census
#     * 2000 block geographies from TIGER/Line
#     * 2010 road geographies from TIGER/Line
#     * NLCD land cover and impervious surface rasters from MRLC
# 2. Rasterize Census blocks
# 3. Delineate "inhabited zones":
#     1. Buffer road geographies by 300 meters and rasterize
#     2. Subtract areas with impervious surface percentage <5% from the
#        rasterized roads
#     3. Subtract water bodies from rasterized roads
# 4. Mask Census blocks to the inhabited zones
# 5. Spread block population evenly across inhabited zone
#
# To combine into a single, nationwide file:
#
# 1. Build a VRT that points to the individual GeoTIFFs
#
#     gdalbuildvrt temp.vrt output/bd-pops/*population_2000.tif
#
# 2. Convert the VRT into a GeoTIFF and clean up
#
#     gdal_translate -of GTiff -co "COMPRESS=LZW" -co "BIGTIFF=YES" temp.vrt output/bd-pops-conus.tif
#     rm temp.vrt
#
# Primary reference material:
#
# 1. IPUMS. (n.d.). 2000 Block Data Standardized to 2010 Geography. IPUMS NHGIS.
# https://www.nhgis.org/documentation/time-series/2000-blocks-to-2010-geog.
# 2. Ruther, M., Leyk, S., & Buttenfield, B. P. (2015). "Comparing the effects
# of an NLCD-derived dasymetric refinement on estimation accuracies for multiple
# areal interpolation methods." GIScience & Remote Sensing 52(2), 158-178.
# http://dx.doi.org/10.1080/15481603.2015.1018856
#
# Contact: Edgar Castro <edgar_castro@g.harvard.edu>

#library(doSNOW)
#library(parallel)
#library(progress)

library(data.table)
library(dplyr)
library(terra)
library(sf)

# Distance, in meters, to buffer around roads in calculating inhabited areas
road_buffer_m <- 300

# Percent impervious area under which land is considered uninhabited
impervious_cutoff <- 5

# Number of cores to use in road buffering
#max_cores <- 20

# Directory to save output to
bd_output_directory <- "output/bd-pops/"

# State FIPS code to produce crosswalk for. Could theoretically do the entire US
# at once but that would cause long periods of unresponsiveness... also maybe
# can be parallelized this way?
#
# This will be overwritten if a state FIPS code is provided as a CLI arg
current_statefp <- "25"

# Very important! We will be generating very large (~1GB) temporary files -
# don't want to use up all of /tmp/ (default location)
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

## State-specific data ----

# 2000 TIGER tracts
tiger_tracts_2000 <- sprintf("/media/qnap3/ShapeFiles/Polygons/TIGER2000/tracts/tl_2010_%s_tract00.shp", current_statefp) %>%
  st_read() %>%
  st_transform(nlcd_crs)

# 2000 TIGER blocks
tiger_blocks_2000 <- sprintf("/media/qnap3/ShapeFiles/Polygons/TIGER2000/blocks/tl_2010_%s_tabblock00.shp", current_statefp) %>%
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
  sprintf("%s/%s_blocks_2000.gpkg", bd_output_directory, current_statefp),
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
  sprintf("%s/%s_blocks_2000.tif", bd_output_directory, current_statefp),
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
  load_function = rast
)

## Buffer and rasterize roads ----
# Here we will be looping over tracts because we can afford to process larger
# areas and move a little faster
#
# UPDATE: see stage1alt_postgis_buffer.R for faster solution using PostGIS

roads_rasterized <- load_cached_data(
  sprintf("%s/%s_roads_buffered.tif", bd_output_directory, current_statefp),
  function() {
    message(sprintf("Buffering roads %s meters", road_buffer_m))
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
  load_function = rast
)

# Parallel version
# roads_rasterized <- load_cached_data(
#   sprintf("output/bd-pops/%s_roads_buffered.tif", current_statefp),
#   function() {
#     message(sprintf("Building SNOW cluster with %d cores", max_cores))
#     cluster <- makeSOCKcluster(max_cores)
#     registerDoSNOW(cluster)
#     on.exit(stopCluster(cluster))
#    
#     bar <- progress_bar$new(
#       sprintf("Buffering roads %s meters :current/:total (:percent) [:bar] eta :eta", road_buffer_m),
#       total = nrow(tiger_roads_2010)
#     )
#     results <- foreach(
#       road = tiger_roads_2010$geometry,
#       .packages = c("progress", "sf"),
#       .options.snow = list(progress = function(n) bar$tick())
#     ) %dopar% {
#       return(st_buffer(road, road_buffer_m))
#     }
# 
#     message("Rasterizing roads")
#     results <- vect(st_as_sfc(results))
#     return(rasterize(
#         results,
#         crop(reference_raster, results)
#     ))
#   },
#   save_function = writeRaster,
#   load_function = rast
# )

## Calculation of inhabited areas ----

inhabited_areas <- load_cached_data(
  sprintf("%s/%s_inhabited_areas_2000.tif", bd_output_directory, current_statefp),
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
  load_function = rast
)

## BD interpolation of population counts ----
  
bd_raster <- load_cached_data(
  sprintf("%s/%s_population_2000.tif", bd_output_directory, current_statefp),
  function() {
    # Start with the inhabited areas from earlier - read from disk so we aren't
    # manipulating the original values
    result <- rast(
      sprintf("%s/%s_inhabited_areas_2000.tif", bd_output_directory, current_statefp)
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
  load_function = rast
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

