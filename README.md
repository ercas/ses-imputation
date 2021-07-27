# ses-imputation

Imputation of missing Census data in years not available in the decennial Census or American Community Survey

## About

This repository contains code used to impute missing years of census data and code used to generate weights and crosswalks between different regionalizations of the same data topic. The latter is necessary because the U.S. Census TIGER/Line boundaries underwent significant changes in 2010 and many of the 2000 geographies do not resemble their 2010 equivalents.

The crosswalk process uses Binary Dasymetric (BD) interpolation combined with a variation of Target-Density Weighting (TDW, Schroeder 2007), building on work from Ruther 2015 and IPUMS.

## References

1. IPUMS. (n.d.). 2000 Block Data Standardized to 2010 Geography. IPUMS NHGIS. https://www.nhgis.org/documentation/time-series/2000-blocks-to-2010-geog. 
2. Ruther, M., Leyk, S., & Buttenfield, B. P. (2015). Comparing the effects of an NLCD-derived dasymetric refinement on estimation accuracies for multiple areal interpolation methods. GIScience & Remote Sensing, 52(2), 158–178. https://doi.org/10.1080/15481603.2015.1018856 
2. Schroeder, J. P. 2007. “Target-Density Weighting Interpolation and Uncertainty Evaluation for Temporal Analysis of Census Data.” Geographical Analysis 39: 311–335. http://dx.doi.org/10.1111/j.1538-4632.2007.00706.x