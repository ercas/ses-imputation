library(dplyr)
library(sf)
library(nabor)

# Load data ---------------------------------------------------------------

tiger_2000_centroids <- st_read("output/zctas/zcta5_2000.gpkg") %>%
  st_centroid()
tiger_2010_centroids <- st_read("output/zctas/zcta5_2010.gpkg") %>%
  st_centroid()
tiger_intersections <- st_read("output/zctas/zcta5_2000_2010_intersections.gpkg") %>%
  st_drop_geometry()

# Find missing geographies ------------------------------------------------

tiger_2000_unmatched <- filter(
  tiger_2000_centroids,
  !ZCTA5CE00 %in% unique(tiger_intersections$ZCTA5CE00)
)
nrow(tiger_2000_unmatched)

tiger_2010_unmatched <- filter(
  tiger_2010_centroids,
  !ZCTA5CE10 %in% unique(tiger_intersections$ZCTA5CE10)
)
nrow(tiger_2010_unmatched)

# Nearest-neighbor match --------------------------------------------------

tiger_2000_matches <- slice(
  tiger_2010_centroids,
  knn(
    st_coordinates(tiger_2010_centroids),
    st_coordinates(tiger_2000_unmatched),
    k = 1
  )$nn.idx[,1]
) %>%
  pull(ZCTA5CE10)

tiger_2010_matches <- slice(
  tiger_2000_centroids,
  knn(
    st_coordinates(tiger_2000_centroids),
    st_coordinates(tiger_2010_unmatched),
    k = 1
  )$nn.idx[,1]
) %>%
  pull(ZCTA5CE00)

rbind(
  data.frame(
    ZCTA5CE00 = tiger_2010_matches,
    ZCTA5CE10 = tiger_2010_unmatched$ZCTA5CE10,
    missing = rep("ZCTA5CE00", nrow(tiger_2010_unmatched)) # Hack to accommodate zero unmatched
  ),
  data.frame(
    ZCTA5CE00 = tiger_2000_matches,
    ZCTA5CE10 = tiger_2000_unmatched$ZCTA5CE00,
    missing = rep("ZCTA5CE10", nrow(tiger_2000_unmatched)) # Hack to accommodate zero unmatched
  )
)  %>%
  write.csv("output/matches.csv", row.names = FALSE)
