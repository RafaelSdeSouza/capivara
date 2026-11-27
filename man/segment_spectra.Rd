\name{segment_spectra}
\alias{segment_spectra}
\title{Segment spectra using a declared scientific representation}
\usage{
segment_spectra(
  cube, n_regions, representation, variance = NULL,
  sample_validity = NULL, support = NULL, redshift = NA_real_,
  wavelength_frame = NULL
)
}
\arguments{
\item{cube}{FITS-like cube or numeric three-dimensional array.}
\item{n_regions}{Positive integer requested hierarchy cut.}
\item{representation}{A non-null representation declaration, such as \code{\link{full_flux_representation}()} or \code{\link{manga_sandra_spectral_shape}()}.}
\item{variance}{Optional aligned variance cube, not inverse variance.}
\item{sample_validity}{Optional aligned logical voxel eligibility array.}
\item{support}{Optional explicit spatial support object.}
\item{redshift}{Systemic redshift, required where rest-frame selection needs it.}
\item{wavelength_frame}{Native wavelength frame, if not declared in cube metadata.}
}
\description{Use the observed-entry spatial Ward hierarchy. Unlike the historical interfaces, this entry point requires an explicit representation and requested region count. Unsupported merges can leave more regions.}
\value{The unchanged CAPIVARA segmentation object returned by \code{\link{segment_large}()}.}
