# elevatr
An R package for accessing elevation data from [Mapzen Tile Service](https://mapzen.com/documentation/terrain-tiles/) and the [Mapzen Elevation Service](https://mapzen.com/documentation/elevation/elevation-service/)

Two services available from Mapzen.com provide access to elevation data as a raster digital elevation model or as height at a single point.  This package provides access to those services and returns elevation data either as a data frame (from the Mapzen Elevation Service) or as a raster object (from the Mapzen Terrain Service).  Additional elevation data sources may be added.

Current plan for this package includes just two functions to access elevation web services:
    - `get_elev_point()`:  Get point elevations using the Mapzen Elevation Service or (for the US Only) the USGS Elevation Point Query Service.  This will accept a data frame of x (long) and y (lat) or a SpatialPoints/SpatialPointsDataFame as input.  A SpatialPointsDataFrame is returned.
    - `get_elev_raster()`: Get elevation data as a raster (e.g. a Digital Elevation Model) from the Mapzen Terrain GeoTIFF Service.  Other sources may be added later.  This will accepts a data frame of of x (long) and y (lat) or any `sp` or `raster` object as input and will return a `raster` object of the elevation tiles that cover the bounding box of the input spatial data. 

## Installation

This package is currently in development and should not be considered stable.  The functions and API may change drastically and rapidly and it may not work at any given moment...  That being said, install with `devtools`


```r
library(devtools)
install_github("jhollist/elevatr")
```

## Attribution
Mapzen terrain tiles contain 3DEP, SRTM, and GMTED2010 content courtesy of the U.S. Geological Survey and ETOPO1 content courtesy of U.S. National Oceanic and Atmospheric Administration.
