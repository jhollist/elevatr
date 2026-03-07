#' function to convert lat long to xyz tile with decimals
#' rounding to tile occurs in \code{get_tilexy}
#' @keywords internal
latlong_to_tilexy <- function(lon_deg, lat_deg, zoom){
  # Code assumes lon is 180 to 180, so converts to that
  lon_deg <- ifelse(lon_deg > 180, lon_deg - 360, lon_deg)
  #Code from https://wiki.openstreetmap.org/wiki/Slippy_map_tilenames#Coordinates_to_tile_numbers_2
  lat_rad <- lat_deg * pi /180
  n <- 2.0 ^ zoom
  xtile <- (lon_deg + 180.0) / 360.0 * n
  ytile <- (1.0 - log(tan(lat_rad) + (1 / cos(lat_rad))) / pi) / 2.0 * n
  return(c(xtile, ytile))
}

#' function to get a data.frame of all xyz tiles to download
#' @keywords internal
get_tilexy <- function(bbx,z){
  #Convert to -180 - +180
  bbx[c("xmin","xmax")] <- ifelse( bbx[c("xmin","xmax")] > 180,  bbx[c("xmin","xmax")] - 360,  bbx[c("xmin","xmax")])
  min_tile <- unlist(slippymath::lonlat_to_tilenum(bbx["xmin"],bbx["ymin"],z))
  max_tile <- unlist(slippymath::lonlat_to_tilenum(bbx["xmax"],bbx["ymax"],z))
  x_all <- seq(from = floor(min_tile[1]), to = floor(max_tile[1]))
  y_all <- seq(from = floor(min_tile[2]), to = floor(max_tile[2]))
  
  
  if(z == 1){
    x_all <- x_all[x_all<2]
    y_all <- y_all[y_all<2]
  } else if(z == 0){
    x_all <- x_all[x_all<1]
    y_all <- y_all[y_all<1]
  }
  
  
  return(expand.grid(x_all,y_all))
}


#' function to check input type and projection.  All input types convert to a
#' SpatialPointsDataFrame for point elevation and bbx for raster.
#' @keywords internal
loc_check <- function(locations, prj = NULL){
 
  
  if(all(class(locations)=="data.frame")){ 
    if(is.null(prj)) { # & !any(class(locations) %in% c("sf", "sfc", "sfg"))){
      stop("Please supply a valid sf crs via argument 'prj'.")
    }    
    locations <- sf::st_as_sf(x = locations, coords = c("x", "y"), crs = prj)
    
  } else if(any(class(locations) %in% c("sf", "sfc", "sfg"))){
    
    sf_crs <- sf::st_crs(locations)
    
    if((is.null(sf_crs) | is.na(sf_crs)) & is.null(prj)){
      stop("Please supply an sf object with a valid crs, or supply a crs with argument 'prj'.")
    }
    
  } else if(any(class(locations) %in% c("SpatRaster", "SpatVector"))){
    
    locations <- sf::st_as_sf(terra::as.points(locations, values = FALSE))
    sf_crs <- sf::st_crs(locations)
    if((is.null(sf_crs) | is.na(sf_crs)) & is.null(prj)){
      stop("Please supply a SpatRaster/Vector with a valid crs, or supply a crs with argument 'prj'.")
    }
  } else if(inherits(locations, "SpatExtent")) {
	if (is.null(prj)) {
		# you could warn, but that seems unnecessary, lonlat is the only reasonable crs in this case 
		# warn("assuming extent has longitude/latitude crs")
		prj <- sf::st_crs(4326)
	}
    locations <- sf::st_as_sf(terra::as.points(locations))
	sf::st_crs(locations) <- prj
  }
  locations$elevation <- 0 

  #check for long>180
  if(is.null(prj)){
    prj_test <- sf::st_crs(locations)
  } else {
    prj_test <- prj
  }
  

  lll <- sf::st_is_longlat(prj_test)
  
  if(lll){
    if(any(sf::st_coordinates(locations)[,1]>180)){
      stop("The elevatr package requires longitude in a range from -180 to 180.")
    } 
  }
  
locations
} 



#' function to project bounding box and if needed expand it
#' @keywords internal
proj_expand <- function(locations,prj,expand){
  
  lll <- sf::st_is_longlat(prj)
    #any(grepl("\\bGEOGCRS\\b",sf::st_crs(prj)) |
    #           grepl("\\bGEODCRS\\b", sf::st_crs(prj)) |
    #           grepl("\\bGEODETICCRS\\b", sf::st_crs(prj)) |
    #           grepl("\\bGEOGRAPHICCRS\\b", sf::st_crs(prj)) |
    #           grepl("\\blonglat\\b", sf::st_crs(prj)) |
    #           grepl("\\blatlong\\b", sf::st_crs(prj)))
  
  if(is.null(nrow(locations))){
    nfeature <- length(locations) 
  } else {
    nfeature <- nrow(locations)
  }
  
  if(any(sf::st_bbox(locations)[c("ymin","ymax")] == 0) & lll & is.null(expand)){
    # Edge case for lat exactly at the equator - was returning NA
    expand <- 0.01
  } else if(nfeature == 1 & lll & is.null(expand)){
    # Edge case for single point and lat long
    expand <- 0.01
  } else if(nfeature == 1 & is.null(expand)){
    # Edge case for single point and projected
    # set to 1000 meters
    unit <- sf::st_crs(sf::st_as_sf(locations), parameters = TRUE)$ud_unit
    expand <- units::set_units(units::set_units(1000, "m"), unit, 
                               mode = "standard")
    expand <- as.numeric(expand)
  }

  
  if(!is.null(expand)){
    
    bbx <- sf::st_bbox(locations) + c(-expand, -expand, expand, expand)
  } else {
    bbx <- sf::st_bbox(locations)
  }

  bbx <- bbox_to_sf(bbx, prj = prj)
  bbx <- sf::st_bbox(sf::st_transform(bbx, crs = ll_geo))
  bbx_coord_check <- as.numeric(bbx)
  if(any(!bbx_coord_check >= -180 & bbx_coord_check <= 360)){
    stop("The elevatr package requires longitude in a range from -180 to 180.")
  } 

  bbx
  
  #sf expand - save for later
  #loc_sf <- sf::st_as_sf(locations)
  #loc_bbx <- sf::st_bbox(loc_sf)
  #bbx_sf <- loc_bbx + c(-expand,-expand,expand,expand)
  #names(bbx_sf) <- c("xmin", "ymin", "xmax", "ymax")
  #attr(bbx_sf, "class") <- "bbox"
  #bbx_sf
}

#' function to clip the DEM
#' @keywords internal
clip_it <- function(rast, loc, expand, clip){
  loc_wm <- sf::st_transform(loc, crs = terra::crs(rast))
  if(clip == "locations" & !grepl("sfc_POINT", class(sf::st_geometry(loc_wm))[1])){
    dem <- terra::mask(terra::crop(rast,loc_wm), loc_wm)
  } else if(clip == "bbox" | grepl("sfc_POINT", class(sf::st_geometry(loc_wm))[1])){
    bbx <- proj_expand(loc_wm, as.character(terra::crs(rast)), expand)
    bbx_sf <- sf::st_transform(bbox_to_sf(bbx), crs = terra::crs(rast))
    dem <- terra::mask(terra::crop(rast,bbx_sf), bbx_sf)
  }
  dem
}

#' Assumes geographic projection
#' sf bbox to poly
#' @param bbox an sf bbox object
#' @param prj defaults to "EPSG:4326"
#' @keywords internal
bbox_to_sf <- function(bbox, prj = 4326) {
  sf_bbx <- sf::st_as_sf(sf::st_as_sfc(bbox))
  sf_bbx <- sf::st_transform(sf_bbx, crs = prj)
  sf_bbx
}

#' Estimate download size of DEMs
#' @param locations the locations
#' @param prj prj string as set earlier by get_elev_point or get_elev_raster
#' @param src the src
#' @param z zoom level if source is aws
#' @keywords internal
estimate_raster_size <- function(locations, prj, src, z = NULL){
  
  locations <- bbox_to_sf(sf::st_bbox(locations), 
                          prj = prj)

  locations <- sf::st_transform(locations, crs = 4326)
  # Estimated cell size (at equator) from zoom level source
  # https://github.com/tilezen/joerd/blob/master/docs/data-sources.md#sources-native-resolution
  # Each degree at equator = 111319.9 meters
  # Convert ground res to dd
  # zoom level 0 = 156543 meters 156543/111319.9
  # old resolution (no idea how I calculated these...)
  # c(0.54905236, 0.27452618, 0.15455633, 0.07145545, 0.03719130, 0.01901903, 
  # 0.00962056, 0.00483847, 0.00241219, 0.00120434, 0.00060173, 0.00030075, 
  #  0.00015035, 0.00007517, 0.00003758)
  m_at_equator <- c(156543.0, 78271.5, 39135.8, 19567.9, 9783.9, 4892.0, 2446.0, 
                    1223.0, 611.5, 305.7, 152.9, 76.4, 38.2, 19.1, 9.6, 4.8, 
                    2.4)
  z_res <- data.frame(z = 0:16, res_dd = m_at_equator/111319.9)
 
  bits <- switch(src,
                 aws = 32,
                 gl3 = 32,
                 gl1 = 32,
                 alos = 32,
                 srtm15plus = 32)
  if(src == "aws"){
    res <- z_res[z_res$z == z,]$res_dd
  } else{
    res <- switch(src,
                  gl3 = 0.0008333,
                  gl1 = 0.0002778,
                  alos = 0.0002778,
                  srtm15plus = 0.004165) 
  }
  num_rows <- (sf::st_bbox(locations)$xmax - sf::st_bbox(locations)$xmin)/res
  num_cols <- (sf::st_bbox(locations)$ymax - sf::st_bbox(locations)$ymin)/res
  
  num_megabytes <- (num_rows * num_cols * bits)/8388608
  num_megabytes
}

#' OpenTopo Key
#' 
#' The OpenTopography API now requires an API Key.  This function will grab your
#' key from an .Renviron file
#' 
#' @keywords internal
get_opentopo_key <- function(){
  if(Sys.getenv("OPENTOPO_KEY")==""){
    stop("You have not set your OpenTopography API Key.
         Please use elevatr::set_opentopo_key().")
  }
  Sys.getenv("OPENTOPO_KEY")
}