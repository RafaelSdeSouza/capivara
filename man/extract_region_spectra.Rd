\name{extract_region_spectra}
\alias{extract_region_spectra}
\alias{plot_regions}
\alias{plot_region_spectra}
\title{Extract and visualize regional spectral products}
\usage{
extract_region_spectra(segmentation, variance = NULL, variance_inflation = 1)
plot_regions(segmentation, ...)
plot_region_spectra(segmentation, ...)
}
\arguments{
\item{segmentation}{A CAPIVARA segmentation result.}
\item{variance}{Optional aligned variance cube.}
\item{variance_inflation}{Positive variance multiplier, default one.}
\item{...}{Arguments to \code{\link{plot_cluster}} or \code{\link{plot_cluster_spectra}}, respectively.}
}
\description{These convenience interfaces retain the existing algorithms, return structures and defaults. Spectral summaries are sums over finite contributors: inspect \code{finite_counts} before treating them as fixed-aperture spectra. Variance propagation is diagonal, with an optional scalar inflation, not a full covariance calculation.}
\details{The region map retains the underlying palette and coordinate defaults. The spectral profile plot retains the historical median-centred display and requires compatible FITS axis metadata; it does not plot calibrated summed flux. Use the extracted spectra to display calibrated sums.}
\value{For extraction, the original spectral summary list; for plotting, the original ggplot object.}
