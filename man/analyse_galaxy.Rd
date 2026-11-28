\name{analyse_galaxy}
\alias{analyse_galaxy}
\title{Analyse a prepared galaxy cube}
\usage{
analyse_galaxy(prepared, n_regions, representation, output_dir = NULL)
}
\arguments{
\item{prepared}{Prepared input list with flux, variance, valid, host_mask, wavelength, wavelength_frame, wavelength_medium, flux_unit, redshift, object_id and provenance. See the Analyse your own data vignette.}
\item{n_regions}{Positive integer requested hierarchy cut.}
\item{representation}{Explicit scientific representation declaration.}
\item{output_dir}{Optional new directory for analysis products. If omitted, the result stays in memory and no files are written.}
}
\description{Define spectral regions and retain raw, quality-filtered and fixed-aperture spectral products. This workflow does not perform physical fitting.}
\value{A capivara_galaxy_analysis object with segmentation, regional spectra, measurements and provenance. Incomplete fixed-aperture channels are missing; uncertainties are diagonal only.}
