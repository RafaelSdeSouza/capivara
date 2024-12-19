#' Plot Median Spectrum and Associated Spectra for a Bin
#'
#' This function visualizes the median spectrum and all associated spectra for a given bin
#' from the output of the Capivara clustering process. It overlays the median spectrum
#' on the plot of all individual spectra for the bin.
#'
#' @param cluster_result A list containing the results from the Capivara clustering process.
#'   Must include the following elements:
#'   \itemize{
#'     \item \code{cluster_map}: A 2D matrix mapping each pixel to a cluster.
#'     \item \code{original_cube}: The original 3D data cube (e.g., IFU data).
#'     \item \code{axDat}: Metadata for the data cube axes, including wavelength information.
#'   }
#' @param bin_id An integer specifying the ID of the bin to be visualized.
#' @param highlight_color A character string specifying the color for the median spectrum
#'   line in the plot (default: \code{"#FFDE38"}).
#'
#' @details
#' The function extracts the spectra for all pixels in the specified bin, computes the median spectrum,
#' and creates a plot where:
#' \itemize{
#'   \item Gray lines represent individual spectra for all pixels in the bin.
#'   \item The highlighted line represents the median spectrum.
#' }
#' The plot is returned as a \code{ggplot2} object for further customization or saving.
#'
#' @return A \code{ggplot2} object representing the plot.
#'
#' @examples
#' # Example usage
#' cluster_result <- list(
#'   cluster_map = matrix(sample(1:3, 100, replace = TRUE), nrow = 10),
#'   original_cube = list(imDat = array(runif(1000), dim = c(10, 10, 10))),
#'   axDat = list(c(1:10, 1:10, seq(4000, 7000, length.out = 10)))
#' )
#'
#' plot <- capivara_plot_spectra(cluster_result, bin_id = 1)
#' print(plot)
#'
#' @import ggplot2
#' @importFrom tidyr pivot_longer
#' @importFrom FITSio axVec
#'
#' @export
capivara_plot_spectra_with_map <- function(cluster_result, bin_id, highlight_color = "#FFDE38", map_palette = "magma") {
  library(ggplot2)
  library(tidyr)
  library(reshape2)
  library(dplyr)
  library(viridis)

  # Extract data from Capivara outputs
  cluster_map <- cluster_result$cluster_map
  original_cube <- cluster_result$original_cube
  wavelengths <- FITSio::axVec(3, cluster_result$axDat)  # Extract wavelength axis

  # Validate bin_id
  if (!bin_id %in% unique(as.vector(cluster_map))) {
    stop(paste("Bin ID", bin_id, "not found in the cluster map."))
  }

  # Extract spectra for the specified bin
  bin_indices <- which(cluster_map == bin_id, arr.ind = TRUE)  # Get pixel coordinates for the bin
  spectra_list <- apply(bin_indices, 1, function(idx) {
    original_cube$imDat[idx[1], idx[2], ]  # Extract spectrum for each pixel in the bin
  })

  # Combine all spectra into a matrix
  spectra_matrix <- t(spectra_list)

  # Compute the median spectrum for the bin
  median_spectrum <- apply(spectra_matrix, 2, median, na.rm = TRUE)
  median_df <- data.frame(
    Wavelength = wavelengths,
    Median_Spectrum = median_spectrum
  )

  # Prepare data for plotting spectra
  spectra_df <- as.data.frame(spectra_matrix)
  spectra_df$Spectrum_ID <- 1:nrow(spectra_matrix)
  spectra_long <- tidyr::pivot_longer(
    spectra_df,
    cols = -Spectrum_ID,
    names_to = "Wavelength_Index",
    values_to = "Flux"
  )
  spectra_long$Wavelength_Index <- as.numeric(gsub("V", "", spectra_long$Wavelength_Index))
  spectra_long$Wavelength <- wavelengths[spectra_long$Wavelength_Index]

  # Prepare data for 2D cluster map
  cluster_df <- reshape2::melt(cluster_map, varnames = c("Row", "Col"), value.name = "Cluster")
  cluster_df <- cluster_df %>%
    mutate(Highlighted = ifelse(Cluster == bin_id, bin_id, NA))  # Highlight the selected bin

  # Plot 1: Spectra
  spectra_plot <- ggplot() +
    geom_line(data = spectra_long, aes(x = Wavelength, y = Flux, group = Spectrum_ID), color = "gray70", size = 0.4) +
    geom_line(data = median_df, aes(x = Wavelength, y = Median_Spectrum), color = highlight_color, size = 1.2) +
    labs(
      title = paste("Spectra for Bin", bin_id),
      x = "Wavelength (Å)",
      y = "Flux"
    ) +
    theme_minimal() +
    theme(panel.grid.minor = element_blank())

  # Plot 2: 2D Cluster Map
  cluster_plot <- ggplot(cluster_df, aes(x = Col, y = Row, fill = as.factor(Highlighted))) +
    geom_tile(color = NA) +
    scale_fill_manual(values = c(NA, highlight_color), na.value = "white") +
    coord_fixed() +
    labs(
      title = paste("2D Map for Bin", bin_id),
      fill = "Highlighted Bin"
    ) +
    theme_void() +
    theme(
      legend.position = "right",
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank()
    )

  # Return both plots as a list
  return(list(spectra_plot = spectra_plot, cluster_plot = cluster_plot))
}
