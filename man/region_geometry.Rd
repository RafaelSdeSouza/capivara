\name{region_geometry}
\alias{region_geometry}
\alias{list_emission_lines}
\title{Inspect region geometry and the optical line registry}
\usage{
region_geometry(segmentation, coords = c("auto", "native", "sky", "physical"))
list_emission_lines(wavelength_medium = c("vacuum", "air"))
}
\arguments{
\item{segmentation}{A CAPIVARA segmentation result.}
\item{coords}{Coordinate convention passed to \code{\link{regions_sf}}.}
\item{wavelength_medium}{Laboratory wavelength convention.}
}
\description{Inspect region geometry and the optical line registry.}
\value{An sf geometry object for region_geometry, or a line registry data frame for list_emission_lines.}
