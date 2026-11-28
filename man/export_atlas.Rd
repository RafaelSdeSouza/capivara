\name{export_atlas}
\alias{export_atlas}
\title{Export an interactive semantic atlas}
\usage{
export_atlas(analysis, output_dir = NULL)
}
\arguments{
\item{analysis}{A capivara_galaxy_analysis object or a path to an existing analysis directory containing the tutorial products.}
\item{output_dir}{New output directory for an in-memory analysis. Omit for an existing analysis directory. Existing files are not overwritten.}
}
\description{Accept an in-memory result from \code{\link{analyse_galaxy}()}, or an existing prepared cube tutorial output directory. No fitting engine is needed. Requires the optional jsonlite package.}
\value{Invisibly, paths to the atlas JSON and two catalogue CSVs.}
