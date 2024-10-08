% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_elev_raster.R
\name{set_opentopo_key}
\alias{set_opentopo_key}
\title{Store OpenTopography Key}
\usage{
set_opentopo_key(key)
}
\arguments{
\item{key}{An OpenTopography API Key as a character.  For details on obtaining an 
OpenTopgraphy key see \url{https://opentopography.org/blog/introducing-api-keys-access-opentopography-global-datasets}.}
}
\description{
This function stores an OpenTopgrapy key in a local .Renviron file. If the
.Renviron file exists, the key will be appended. This will typically only 
need to be done once per machine.
}
