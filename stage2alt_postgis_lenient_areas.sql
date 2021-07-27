-- PostGIS code to calculate the area within 300 meters of roads intersecting
-- each ZCTA. Assumes that the SQL code in stage1alt_rasterize_postgis_buffer.R
-- has already been run.
-- $ shp2pgsql -I -s 4269 tl_2010_us_zcta510.shp tl_2010_us_zcta510 | psql postgis
--
-- Contact: Edgar Castro <edgar_castro@g.harvard.edu>

-- Intersection of road buffers with 2000 ZCTA5 geographies
CREATE TABLE tl_2010_us_roads_zcta500_intersections AS
    SELECT tl_2010_us_zcta500_valid.zcta5ce00,
           tl_2010_us_roads_300mbuffer.gid,
           ST_Multi(ST_Buffer(ST_Intersection(tl_2010_us_zcta500_valid.geom, tl_2010_us_roads_300mbuffer.geom), 0.0)) AS geom
    FROM tl_2010_us_zcta500_valid
    INNER JOIN tl_2010_us_roads_300mbuffer
        ON ST_Intersects(tl_2010_us_zcta500_valid.geom, tl_2010_us_roads_300mbuffer.geom)
    WHERE NOT ST_IsEmpty(ST_Buffer(ST_Intersection(tl_2010_us_zcta500_valid.geom, tl_2010_us_roads_300mbuffer.geom), 0.0));
    
-- Merged version
CREATE TABLE tl_2010_us_roads_zcta500_intersections_merged AS
    SELECT zcta5ce00, ST_Union(geom)
    FROM tl_2010_us_roads_zcta500_intersections
    GROUP BY zcta5ce00;
    
-- Intersection of road buffers with 2010 ZCTA5 geographies
CREATE TABLE tl_2010_us_roads_zcta510_intersections AS
    SELECT tl_2010_us_zcta510_valid.zcta5ce10,
           tl_2010_us_roads_300mbuffer.gid,
           ST_Multi(ST_Buffer(ST_Intersection(tl_2010_us_zcta510_valid.geom, tl_2010_us_roads_300mbuffer.geom), 0.0)) AS geom
    FROM tl_2010_us_zcta510_valid
    INNER JOIN tl_2010_us_roads_300mbuffer
        ON ST_Intersects(tl_2010_us_zcta510_valid.geom, tl_2010_us_roads_300mbuffer.geom)
    WHERE NOT ST_IsEmpty(ST_Buffer(ST_Intersection(tl_2010_us_zcta510_valid.geom, tl_2010_us_roads_300mbuffer.geom), 0.0));
    
-- Merged version
CREATE TABLE tl_2010_us_roads_zcta510_intersections_merged AS
    SELECT zcta5ce10, ST_Union(geom)
    FROM tl_2010_us_roads_zcta510_intersections
    GROUP BY zcta5ce10;

-- Intersection of road buffers with 2000-2010 ZCTA5 intersection geographies
CREATE TABLE tl_2010_us_roads_zcta500_zcta510_intersections_intersections AS
    SELECT tl_2010_us_zcta500_zcta510_intersections.zcta5ce10,
           tl_2010_us_zcta500_zcta510_intersections.zcta5ce00,
           tl_2010_us_roads_300mbuffer.gid,
           ST_Multi(ST_Buffer(ST_Intersection(tl_2010_us_zcta500_zcta510_intersections.geom, tl_2010_us_roads_300mbuffer.geom), 0.0)) AS geom
    FROM tl_2010_us_zcta500_zcta510_intersections
    INNER JOIN tl_2010_us_roads_300mbuffer
        ON ST_Intersects(tl_2010_us_zcta500_zcta510_intersections.geom, tl_2010_us_roads_300mbuffer.geom)
    WHERE NOT ST_IsEmpty(ST_Buffer(ST_Intersection(tl_2010_us_zcta500_zcta510_intersections.geom, tl_2010_us_roads_300mbuffer.geom), 0.0));

-- Merged version
CREATE TABLE tl_2010_us_roads_zcta500_zcta510_intersections_intersections_merged AS
    SELECT zcta5ce10, zcta5ce00, ST_Union(geom)
    FROM tl_2010_us_roads_zcta500_zcta510_intersections_intersections
    GROUP BY zcta5ce10, zcta5ce00;

-- Export merged versions
-- $ pgsql2shp -f output/zctas/tl_2010_us_roads_zcta500_intersections_merged.shp -h localhost postgis tl_2010_us_roads_zcta500_intersections_merged
-- $ pgsql2shp -f output/zctas/tl_2010_us_roads_zcta510_intersections_merged.shp -h localhost postgis tl_2010_us_roads_zcta510_intersections_merged
-- $ pgsql2shp -f output/zctas/tl_2010_us_roads_zcta500_zcta510_intersections_intersections_merged.shp -h localhost postgis tl_2010_us_roads_zcta500_zcta510_intersections_intersections_me