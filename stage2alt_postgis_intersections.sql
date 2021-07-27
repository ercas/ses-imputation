-- PostGIS code to generate ZCTA5 2000 <-> 2010 intersections - much faster than
-- an sf solution due to the persistent spatial index
--
-- Timing information:
-- 0.00user 0.00system 7:07.31elapsed 0%CPU (0avgtext+0avgdata 7348maxresident)k
-- 0inputs+0outputs (0major+400minor)pagefaults 0swaps
--
-- Should first import both the zcta500 and zcta510 shapefiles from TIGER2010:
-- $ shp2pgsql -I -s 4269 tl_2010_us_zcta500.shp tl_2010_us_zcta500 | psql postgis
-- $ shp2pgsql -I -s 4269 tl_2010_us_zcta510.shp tl_2010_us_zcta510 | psql postgis
--
-- Contact: Edgar Castro <edgar_castro@g.harvard.edu>

-- Validate, index, and vacuum analyze ZCTA5 2000 geographies
CREATE TABLE tl_2010_us_zcta500_valid AS
    SELECT zcta5ce00, ST_MakeValid(geom) AS geom
    FROM tl_2010_us_zcta500;
CREATE INDEX tl_2010_us_zcta500_valid_geom_idx
    ON tl_2010_us_zcta500_valid USING GIST(geom);
CREATE INDEX tl_2010_us_zcta500_valid_zcta5ce00_idx
    ON tl_2010_us_zcta500_valid(zcta5ce00);
VACUUM ANALYZE tl_2010_us_zcta500_valid;

-- Validate, index, and vacuum analyze ZCTA5 2010 geographies
CREATE TABLE tl_2010_us_zcta510_valid AS
    SELECT zcta5ce10, ST_MakeValid(geom) AS geom
    FROM tl_2010_us_zcta510;
CREATE INDEX tl_2010_us_zcta510_valid_geom_idx
    ON tl_2010_us_zcta510_valid USING GIST(geom);
CREATE INDEX tl_2010_us_zcta500_valid_zcta5ce10_idx
    ON tl_2010_us_zcta510_valid(zcta5ce10);
VACUUM ANALYZE tl_2010_us_zcta510_valid;

-- Create intersections
-- See https://postgis.net/docs/ST_Intersection.html, Examples section for more
-- information about ST_Buffer(geometry, 0.0)
CREATE TABLE tl_2010_us_zcta500_zcta510_intersections AS
    SELECT tl_2010_us_zcta510_valid.zcta5ce10,
           tl_2010_us_zcta500_valid.zcta5ce00,
           ST_Multi(ST_Buffer(ST_Intersection(tl_2010_us_zcta510_valid.geom, tl_2010_us_zcta500_valid.geom), 0.0)) AS geom
    FROM tl_2010_us_zcta510_valid
    INNER JOIN tl_2010_us_zcta500_valid
        ON ST_Intersects(tl_2010_us_zcta510_valid.geom, tl_2010_us_zcta500_valid.geom)
    WHERE NOT ST_IsEmpty(ST_Buffer(ST_Intersection(tl_2010_us_zcta510_valid.geom, tl_2010_us_zcta500_valid.geom), 0.0));
DELETE FROM

-- Now, we can export the new intersections table to a shapefile:
-- $ pgsql2shp -f output/zctas/zcta5_2000_2010_intersections.shp -h localhost postgis tl_2010_us_zcta500_zcta510_intersections