#' Plot a Cluster Map Color-Coded by Signal-to-Noise Ratio (SNR)
#'
#' This function visualizes the output of a clustering algorithm applied to
#' an IFU data cube, where the clusters are color-coded by their respective
#' signal-to-noise ratio (SNR).
#'
#' @param cluster_data A list containing the clustering results with the following components:
#'   \itemize{
#'     \item \code{cluster_map}: A matrix representing the spatial layout of the clusters.
#'     \item \code{cluster_snr}: A numeric vector with the SNR values for each cluster.
#'   }
#' @param palette Character. Name of the Viridis color palette to use for the SNR scale.
#'   Options include "magma", "inferno", "plasma", "viridis", and "cividis". Defaults to "magma".
#'
#' @return A \code{ggplot2} object representing the cluster map, color-coded by SNR.
#'
#' @details
#' The function converts the \code{cluster_map} into a long-format data frame and merges it
#' with the \code{cluster_snr} values. The SNR for each cluster is displayed using a Viridis
#' color scale, ensuring high perceptual uniformity.
#'
#' @examples
#' \dontrun{
#' # Example cluster data
#' cluster_data <- list(
#'   cluster_map = matrix(sample(1:5, 100, replace = TRUE), nrow = 10, ncol = 10),
#'   cluster_snr = c(10, 20, 30, 40, 50)
#' )
#'
#' # Plot the cluster map
#' plot <- plot_cluster(cluster_data, title = "Cluster Map by SNR", palette = "viridis")
#' print(plot)
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_tile coord_fixed labs theme_void theme element_blank
#' @importFrom reshape2 melt
#' @importFrom dplyr mutate
#' @importFrom viridis scale_fill_viridis
#' @export
plot_cluster <- function(cluster_data, palette = "magma") {
  # Extract cluster map and SNR values
  cluster_map <- cluster_data$cluster_map
  cluster_snr <- cluster_data$cluster_snr

  # Convert the cluster map to a long-format data frame
  cluster_df <- reshape2::melt(cluster_map, varnames = c("Row", "Col"), value.name = "Cluster")

  # Optionally add SNR information (here not used for fill)
  cluster_df <- dplyr::mutate(cluster_df, SNR = ifelse(!is.na(Cluster), cluster_snr[Cluster], NA))

  # Convert Cluster to factor for discrete coloring
  cluster_df$Cluster <- as.factor(cluster_df$Cluster)

  # Determine number of clusters (levels)
  n_clusters <- length(levels(cluster_df$Cluster))

  # Create a palette based on the input:
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

  # Generate the plot using a discrete fill scale for clusters
  plot <- ggplot2::ggplot(cluster_df, ggplot2::aes(x = Row, y = Col, fill = Cluster)) +
    ggplot2::geom_tile() +
    ggplot2::coord_fixed() +
    ggplot2::scale_fill_manual(values = colors, na.value = "black") +
    ggplot2::theme_void() +
    ggplot2::theme(
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      legend.position = "none"
    )

  return(plot)
}
