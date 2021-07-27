# Instead of buffering in R, offload the task to PostGIS which can get the job
# done a lot more quickly
#
# 1. Import the NLCD spatial reference system (EPSG:42303, Albers equal area
#    conical, continental U.S., NAD83)
#
#     INSERT into spatial_ref_sys (srid, auth_name, auth_srid, proj4text, srtext) values ( 42303, 'EPSG', 42303, '+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs ', 'PROJCS["NAD83 / Albers NorthAm",GEOGCS["NAD83",DATUM["North_American_Datum_1983",SPHEROID["GRS_1980",6378137,298.257222101]],PRIMEM["Greenwich",0],UNIT["Decimal_Degree",0.0174532925199433]],PROJECTION["Albers_conic_equal_area"],PARAMETER["central_meridian",-96.0],PARAMETER["latitude_of_origin",23],PARAMETER["standard_parallel_1",29.5],PARAMETER["standard_parallel_2",45.5],PARAMETER["false_easting",0],PARAMETER["false_northing",0],UNIT["Meter",1],AUTHORITY["EPSG","42303"]]');
#
# 2. Transform roads to EPSG:42303
#
#     CREATE TABLE tl_2010_us_roads_epsg42303 AS SELECT * FROM tl_2010_us_roads;
#     ALTER TABLE tl_2010_us_roads_epsg42303 ALTER COLUMN geom TYPE Geometry(MultiLineString, 42303) USING ST_Transform(geom, 42303);
#
# 3. Buffer roads
#
#     CREATE TABLE tl_2010_us_roads_epsg42303_300mbuffer AS SELECT * FROM tl_2010_us_roads_epsg42303;
#     ALTER TABLE tl_2010_us_roads_epsg42303_300mbuffer ALTER COLUMN geom TYPE Geometry(MultiPolygon, 42303) USING ST_Multi(ST_Buffer(geom, 300));
#
# 4. Reproject to WGS84 (most programs do not support EPSG:42303)
#
#     CREATE TABLE tl_2010_us_roads_epsg4326_300mbuffer AS SELECT * FROM tl_2010_us_roads_epsg42303_300mbuffer;
#     ALTER TABLE tl_2010_us_roads_epsg4326_300mbuffer ALTER COLUMN geom TYPE Geometry(MultiPolygon, 4326) USING ST_Transform(geom, 4326);
#
# 4.1. Optional: also reproject to 4269 to match the rest of TIGER/Line
#
#     CREATE TABLE tl_2010_us_roads_300mbuffer AS SELECT * FROM tl_2010_us_roads_epsg4326_300mbuffer;
#     ALTER TABLE tl_2010_us_roads_300mbuffer ALTER COLUMN geom TYPE Geometry(MultiPolygon, 4269) USING ST_Transform(geom, 4269);
#
# 5. Export resulting files (Bash)
#
#     for statefp in $(echo {01..99})
#     do
#       echo $statefp
#       pgsql2shp -f ${statefp}_roads_buffered.shp -h 127.0.0.1 postgis \
#         "SELECT * FROM tl_2010_us_roads_epsg4326_300mbuffer WHERE statefp = '$statefp'"
#       if ! [ -f ${statefp}*.shp ]
#       then
#         rm -v ${statefp}*
#       fi
#     done
#
# Contact: Edgar Castro <edgar_castro@g.harvard.edu>

library(sf)
library(stringr)
library(terra)

statefps <- str_extract(Sys.glob("output/roads-buffered/*.shp"), "[0-9]+")

args <- commandArgs(trailingOnly=TRUE)
if (length(args) > 0) {
  if (!args[1] %in% statefps) {
    stop(sprintf("Unsupported state FIPS code %s", args[1]))
  }
  statefps <- args[1]
}

road_buffer_m <- 300

nlcd_land_cover_2001 <- rast("/media/qnap4/MRLC NLCD/landcover/nlcd_2001_land_cover_l48_20210604.img")
nlcd_crs <- crs(nlcd_land_cover_2001)
    
for (current_statefp in statefps) {
  output_file <- sprintf("output/bd-pops/%s_roads_buffered.tif", current_statefp)
  if (file.exists(output_file)) {
    message(sprintf("Skipping: state with FIPS code %s", current_statefp))
  } else {
    message(sprintf("Processing: state with FIPS code %s", current_statefp))
    
    message("> Reading shapefile")
    roads_to_rasterize <- st_read(sprintf("output/roads-buffered/%s_roads_buffered.shp", current_statefp))
    
    message("> Reprojecting to NLCD CRS")
    roads_to_rasterize <- vect(st_transform(roads_to_rasterize, nlcd_crs))
    
    tryCatch(
      {
        message("> Rasterizing buffered roads")
        roads_raster <- rasterize(
            roads_to_rasterize,
            crop(nlcd_land_cover_2001, roads_to_rasterize)
        )
        
        message("> Writing rasterized roads")
        writeRaster(roads_raster, sprintf("output/bd-pops/%s_roads_buffered.tif", current_statefp))
      },
      error = function(x) {
        message(sprintf("State FIPS code %s not  in CONUS", args[1]))
      }
    )
    
  }
}
