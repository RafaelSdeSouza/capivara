\name{explore_galaxy}
\alias{explore_galaxy}
\title{Explore a galaxy in a local interactive atlas}
\usage{
explore_galaxy(analysis, output_dir = tempfile("capivara-atlas-"), open = interactive())
}
\arguments{
\item{analysis}{A result from \code{\link{analyse_galaxy}()}.}
\item{output_dir}{New directory for atlas and reproducibility products.}
\item{open}{Open the HTML in the default browser? Defaults to interactive().}
}
\description{Exports the native region products and embeds the JSON into a local HTML viewer. No data are uploaded. The viewer supports native pixel coordinates and descriptive spectral measurements; it does not infer physical labels.}
\value{Invisibly, the absolute path to the HTML viewer.}
