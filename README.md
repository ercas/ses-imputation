# ses-imputation

Imputation of missing Census data in years not available in the decennial Census or American Community Survey

## About

This repository contains code used to impute missing years of census data and code used to generate weights and crosswalks between different regionalizations of the same data topic. The latter is necessary because the U.S. Census TIGER/Line boundaries underwent significant changes in 2010 and many of the 2000 geographies do not resemble their 2010 equivalents.

The crosswalk process uses Binary Dasymetric (BD) refinemant combined with a variation of Target-Density Weighting (TDW, Schroeder 2007), building on work from Ruther 2015 and IPUMS.

## Process

#### Stage 1: Binary Dasymetric (BD) refinemant of inhabited zones

In the first stage, "inhabited zones" are created from a combination of 2010 U.S. Census Bureau [TIGER/Line shapefiles](https://www.census.gov/geographies/mapping-files/time-series/geo/tiger-line-file.html) and the 2001 MRLC [NLCD percent developed imperviousness](https://www.mrlc.gov/data/nlcd-2001-percent-developed-imperviousness-conus) and [NLCD land cover](https://www.mrlc.gov/data/nlcd-2001-land-cover-conus) rasters, following the methodology from IPUMS.

Inhabited zones are calculated as follows:

1. Buffer the 2010 TIGER/Line road lines by 300 meters;
2. Rasterize these to the 2001 NLCD 30-meter grid (shared by both rasters);
3. Subtract out water bodies from the 2001 NLCD land cover raster;
4. Subtract out grid cells with <5% imperviousness from the 2001 NLCD percent developed imperviousness raster.

### Stage 2: Cascaded Binary Dasymetric Target-Density Weighting (C-BD-TDW) of disparate regionalizations

The second stage generates weights via the TDW method described in Schroeder, 2007, incorporating the inhabited zones from Stage 1 in a manner similar to both IPUMS and Ruther et al., 2015. The original TDW formula is given in Schroeder, 2007 as:

![](static/tdw.svg) (Eq. 1)

where:

* A refers to the area of a geography or intersection;
* Z refers to the population of a geography;
* y refers to the variable being crosswalked;
* Subscript s refers to the "source" universe, i.e. that containing the original data (the 2000 TIGER/Line shapefiles);
* Subscripts t and 𝜏 refer to the "target" universe, i.e. that which the original data is being crosswalked *to* (the 2010 TIGER/Line shapefiles); and
* The combination of subscripts st or s𝜏 refers to the intersection of geographies between the source and target regionalizations.

We slightly depart from Schroeder by using the areas of the inhabited zones produced in Stage 1 as follows:

![](static/tdw2.svg) (Eq. 2)

where A\* now refers to the BD-refined areas rather than the total areas of specific geographies. This can be thought of as Binary-Dasmetric Target-Density Weighting (BD-TDW). In this way, we are able to generate weights that are more representative of the true geospatial distribution of population than would be possible using TIGER/Line alone.

At this point, we now have an issue: some geographies will not have any NLCD grid cells with more than 5% impervious surface, such as those in rural areas. To overcome this, we follow a cascaded BD-TDW (C-BD-TDW) methodology resembling the methodology of IPUMS.

In C-BD-TDW, we first start with the strictest definition of inhabited zones - presented above - in the calculation of A\*. In geographies where A\* is zero due to imperviousness restrictions, we again redefine the TDW formula as

![](static/tdw3.svg) (Eq. 3)

where A<sup>r</sup> now refers to a definition of inhabited zones involving only the 300-meter road buffers, without any restriction on land cover.

In cases where A<sup>r</sup> is zero due to no roads, we finally fall back to the original TDW formula (Eq. 1).

In short, in stage 3 we combine the methodologies of IPUMS and Ruther et al., 2015. in that we use only the raw TDW formula with the BD-refined inhabited zones, as per Ruther et al., 2015, and that we are use the inhabited zones definition and cascading of progressively less-strict inhabited zone definitions as per IPUMS.

To save on processing time, we store the C-BD-TDW outputs separately, consisting of weights (the part in the equations between the summation and y<sub>s</sub>) and GEOIDs of the corresponding geographies in the source and target universes. These can then be reused to translate each geography characteristic separately.

### Stage 3: Matching of remaining geographies

A few TIGER/Line 2000 geographies have zero overlap with 2010 geographies due to significant changes between regionalizations, and vice versa. There is only a small number of geographies with tihs issue, but for completeness, we opt to match them to existing geographies.

In these cases, the 2000 geographies are matched to the 2010 geographies via k-nearest-neighbour with k=1 and straight-line distance on centroids. The same process is repeated to match 2010 geographies to 2000 geographies.

**Possible future improvements**:

* Use a cascaded approach to finding centroids in a manner similar to in C-BD-TDW, using first the centroid of inhabited areas, then the centroid of 300-meter road buffers, then the centroid of the overall geography
* Use a distance function to discard unmatched geographies for which there are no nearby geographies with 

### Stage 4: Interpolation of Census data for missing years (C-BD-TDW)

In the third stage, we can finally interpolate Census data from years not covered by the Decennial Census or American Community Survey.

The first step of the interpolation is ensuring that all data exists in the same regionalization. In the previous two stages, we have facilitated this by generating the weights and crosswalks that translate data from the 2000 TIGER regionalization to the 2010 TIGER regionalization. The translated data can be obtained simply by substituting the weights back into the TDW equation.

If characteristics are in terms of percentages, we then need to do some extra math to account for this, instead using the equation

![](static/tdwpct.svg) (Eq. 4)

where w<sub>s</sub> refers to the weights generated in Stage 2 and Z<sub>s</sub> refers to the population denominator of the variable being translated y<sub>s</sub>, i.e. the total population for population characteristics, the number of homes for housing characteristics, etc. Essentially, this is the TDW formula where the percentage is translated back into a count in the numerator and then translated back into a percentage by translating the population denominator itself.

From here, we then linearly interpolate the missing via simple weighted least-squares regression models where the year is the sole predictor of the geography characteristic, which is the response variable. Data from different Census data sets are weighted as follows:

| Source                                     | Weight |
|--------------------------------------------|--------|
| Decennial Census                           | 4      |
| American Community Survey, 5-year estimate | 2      |
| American Community Survey, 1-year estimate | 1      |

## Diagnostics

TODO

## References

1. IPUMS. (n.d.). 2000 Block Data Standardized to 2010 Geography. IPUMS NHGIS. https://www.nhgis.org/documentation/time-series/2000-blocks-to-2010-geog. 
2. Ruther, M., Leyk, S., & Buttenfield, B. P. (2015). Comparing the effects of an NLCD-derived dasymetric refinement on estimation accuracies for multiple areal interpolation methods. GIScience & Remote Sensing, 52(2), 158–178. https://doi.org/10.1080/15481603.2015.1018856 
2. Schroeder, J. P. 2007. “Target-Density Weighting Interpolation and Uncertainty Evaluation for Temporal Analysis of Census Data.” Geographical Analysis 39: 311–335. http://dx.doi.org/10.1111/j.1538-4632.2007.00706.x
