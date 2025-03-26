#' Visualize Scaled Spectra and Median Profiles for Clustered IFU Data
#'
#' This function creates a faceted plot that displays the individual pixel spectra (after
#' scaling) and the corresponding median spectrum for each segmented cluster obtained from
#' an IFU data cube. Each facet corresponds to a unique cluster, allowing for easy comparison
#' of spectral profiles across different spatial regions.
#'
#' @param cluster_result A list produced by a segmentation function (e.g. \code{\link{segment}})
#'   containing at least the following elements:
#'   \itemize{
#'     \item \code{cluster_map}: A matrix mapping each spatial pixel to a cluster.
#'     \item \code{original_cube}: The original data cube (typically a FITS object).
#'     \item \code{axDat}: Axis metadata (including wavelength information).
#'   }
#' @param scale_fn A function that scales each individual spectrum (a numeric vector).
#'   The default is \code{median_scale} which should be defined elsewhere in your package.
#' @param highlight_color A character string specifying the color for plotting the median spectrum.
#'   Defaults to \code{"red3"}.
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Extracts the spatial cluster map, original data cube, and wavelength information from
#'         \code{cluster_result}.
#'   \item Converts the cube to a 2D matrix (assuming \code{cube_to_matrix} is available).
#'   \item Creates a data frame mapping each valid spatial pixel to its cluster.
#'   \item Extracts each pixel's spectrum and applies the scaling function \code{scale_fn} to it.
#'   \item Computes the median spectrum for each cluster at every wavelength.
#'   \item Constructs a faceted ggplot where:
#'     \itemize{
#'       \item Gray lines represent the scaled individual spectra.
#'       \item A highlighted line (using \code{highlight_color}) represents the median spectrum.
#'     }
#' }
#'
#' @return A \code{ggplot2} object representing a faceted plot of the scaled individual spectra
#'   and their median profiles for each cluster.
#'
#' @seealso \code{\link{segment}}, \code{\link{cube_to_matrix}}, \code{\link[FITSio]{axVec}}
#'
#' @examples
#' \dontrun{
#'   # Assuming 'input_cube' is a FITS cube and 'segment' clusters the data:
#'   cluster_result <- segment(input_cube, Ncomp = 5)
#'
#'   # Generate the spectra plot with scaling and a highlighted median spectrum
#'   spectra_plot <- plot_cluster_spectra(cluster_result, scale_fn = median_scale, highlight_color = "red3")
#'   print(spectra_plot)
#' }
#'
#' @export
plot_cluster_spectra <- function(cluster_result, scale_fn = median_scale, palette = "magma") {
  # Extract relevant data from the segmentation result
  cluster_map <- cluster_result$cluster_map
  original_cube <- cluster_result$original_cube
  wavelengths <- FITSio::axVec(3, cluster_result$axDat)  # Assumes axDat holds wavelength info

  # Convert the cube to a 2D matrix (assuming cube_to_matrix is available)
  IFU2D <- cube_to_matrix(original_cube)

  # Prepare a data frame with spatial indices and corresponding cluster assignments
  pixel_df <- data.frame(
    index = 1:length(as.vector(cluster_map)),
    cluster = as.vector(cluster_map)
  )
  pixel_df <- pixel_df %>% filter(!is.na(cluster))

  # For each valid pixel, extract its spectrum and store with spatial info
  spectra_list <- lapply(1:nrow(pixel_df), function(i) {
    # Convert linear index to matrix row/col; note: R is column-major, so adjust if needed
    pos <- which(cluster_map == pixel_df$cluster[i], arr.ind = TRUE)[1, ]
    # Apply scaling function to the spectrum if needed:
    raw_spectrum <- original_cube$imDat[pos[1], pos[2], ]
    scaled_spectrum <- scale_fn(raw_spectrum)
    data.frame(
      Wavelength = wavelengths,
      Flux = scaled_spectrum,
      pixel = paste(pos[1], pos[2], sep = ","),
      cluster = factor(pixel_df$cluster[i])
    )
  })
  spectra_df <- do.call(rbind, spectra_list)

  # Compute the median spectrum per cluster & wavelength
  median_spectra <- spectra_df %>%
    group_by(cluster, Wavelength) %>%
    summarise(Median_Flux = median(Flux, na.rm = TRUE), .groups = "drop")


  # Force clusters to be ordered numerically
  median_spectra$cluster <- factor(
    median_spectra$cluster,
    levels = sort(as.numeric(unique(median_spectra$cluster)))
  )

  # Determine the number of clusters
  n_clusters <- length(unique(median_spectra$cluster))

  # Create a discrete palette based on the input palette argument
  if (is.character(palette) && length(palette) == 1) {
    # Generate discrete colors using viridis if a single option is given
    colors <- viridis::viridis(n_clusters, option = palette)
  } else if (is.character(palette) && length(palette) > 1) {
    # If a custom color vector is provided, ensure there are enough colors
    if(length(palette) < n_clusters) {
      colors <- grDevices::colorRampPalette(palette)(n_clusters)
    } else {
      colors <- palette[1:n_clusters]
    }
  } else {
    stop("Palette must be either a single viridis option string or a vector of colors")
  }

  # Plot median spectra colored by cluster
  p <- ggplot2::ggplot(median_spectra, ggplot2::aes(x = Wavelength, y = Median_Flux, color = cluster)) +
    geom_line(linewidth = 0.15) +
    ggplot2::scale_color_manual(values = colors) +
    ggplot2::facet_wrap(~ cluster, scales = "free_y") +
    ggplot2::labs(x = "Wavelength (Å)", y = "Flux") +
    ggplot2::theme(
      panel.background = element_rect(fill = "gray50"),
      axis.ticks.length = unit(0.1, "cm"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "none",
      legend.background = element_rect(colour = "transparent", fill = "transparent"),
      legend.text = element_text(size = 10),
      legend.spacing.x = unit(0.1, 'cm'),
      legend.key = element_rect(colour = "transparent", fill = "transparent"),
      axis.title = element_text(color = "black", size = 11)
    )

  return(p)
}

