#' Approximate Ward Segmentation of a Data Cube (Scalable)
#'
#' This function performs a scalable, approximate version of Ward hierarchical
#' clustering for IFU/hyperspectral data cubes. The cube is flattened into a
#' 2D matrix (spatial pixels \eqn{\times} spectral variables), robustly scaled
#' row-wise, optionally projected into a low-dimensional PCA space, and then
#' clustered using Ward's D2 linkage on a random subset of pixels. All remaining
#' valid pixels are assigned to the nearest cluster centroid in PCA space.
#'
#' Compared to \code{\link{segment}} (exact Ward on all valid pixels), this method
#' avoids constructing an \eqn{O(n^2)} distance object and is therefore suitable
#' for large cubes with many spatial pixels.
#'
#' @param input A FITS object representing the input data cube. Typically, this is an IFU data cube.
#' @param Ncomp Integer, the number of clusters to form.
#' @param sample_n Integer, number of valid pixels sampled to build the Ward dendrogram.
#'   Larger values better approximate exact Ward but increase runtime and memory (\eqn{O(sample\_n^2)}).
#' @param pca_k Integer, number of principal components retained for clustering and assignment.
#'   Must be \eqn{\le} the number of spectral variables. Smaller values can improve robustness and speed,
#'   but may smooth over narrow spectral features.
#' @param redshift Numeric, the redshift to apply for wavelength correction. Defaults to 0 (no correction).
#'
#' @details
#' Steps performed by the function:
#' \enumerate{
#'   \item Reads the input FITS data cube.
#'   \item Converts the cube into a 2D matrix (spatial pixels x spectral variables).
#'   \item Identifies valid pixels by requiring: (i) a minimum fraction of finite spectral values,
#'         (ii) non-zero robust scatter (MAD) across the spectrum, and (iii) positive finite energy.
#'   \item Scales each valid spectrum robustly (median/MAD; falls back to SD if needed), replacing any
#'         non-finite scaled values with zero.
#'   \item Computes a truncated PCA of the scaled spectra using \code{\link[irlba]{prcomp_irlba}} and
#'         retains \code{pca_k} components.
#'   \item Randomly samples \code{sample_n} valid pixels in PCA space, computes pairwise distances,
#'         and performs Ward's D2 hierarchical clustering via \code{\link[fastcluster]{hclust}}.
#'   \item Cuts the dendrogram into \code{Ncomp} clusters and computes cluster centroids in PCA space
#'         from the sampled pixels.
#'   \item Assigns all valid pixels to the nearest centroid in PCA space (blockwise for memory efficiency).
#'   \item Reshapes assignments into a 2D cluster map matching the spatial dimensions of the input cube.
#' }
#'
#' @return A list containing:
#' \item{cluster_map}{A \code{n_rows x n_cols} matrix of cluster assignments (integers), aligned with the original cube layout.}
#' \item{valid}{Integer vector of linear indices for pixels used in the segmentation (after validity filtering).}
#' \item{pca}{The PCA object returned by \code{\link[irlba]{prcomp_irlba}} (useful for diagnostics and projections).}
#' \item{centers}{A \code{Ncomp x pca_k} matrix of cluster centroids in PCA space.}
#'
#' @section Notes on approximation:
#' This function is not identical to exact Ward clustering on all pixels. The dendrogram
#' is built on a subset of pixels and the final labeling is obtained by nearest-centroid
#' assignment in PCA space. For large cubes, this trade-off typically yields similar
#' large-scale segmentations at a fraction of the computational cost.
#'
#' @seealso \code{\link{segment}}, \code{\link[irlba]{prcomp_irlba}},
#'   \code{\link[fastcluster]{hclust}}, \code{\link[stats]{cutree}}
#'
#' @examples
#' \dontrun{
#' # Read a FITS cube and run scalable Ward segmentation
#' input_cube <- FITSio::readFITS("manga_7443_12703_LOGCUBE.fits")
#'
#' seg <- segment_ward_approx(
#'   input    = input_cube,
#'   Ncomp    = 5,
#'   sample_n = 10000,
#'   pca_k    = 20
#' )
#'
#' cluster_map <- seg$cluster_map
#' valid_idx   <- seg$valid
#' }
#'
#' @export
segment_approx  <- function(input, Ncomp = 5, sample_n = 10000,
                                       pca_k = 20, redshift = 0) {
  cubedat <- input
  n_row <- dim(cubedat$imDat)[1]
  n_col <- dim(cubedat$imDat)[2]
  n_wave <- dim(cubedat$imDat)[3]
  
  IFU2D <- cube_to_matrix(cubedat)  # (n_pix x n_wave)
  n_pix <- nrow(IFU2D)
  
  # ---- valid pixels (similar to yours)
  n_finite <- rowSums(is.finite(IFU2D))
  finite_ok <- n_finite >= max(10L, floor(0.8 * n_wave))
  
  row_scatter <- apply(IFU2D, 1, function(v) {
    vv <- v[is.finite(v)]
    if (length(vv) < 2) return(0)
    stats::mad(vv, center = stats::median(vv), constant = 1)
  })
  scatter_ok <- row_scatter > 0
  
  row_energy <- rowSums(IFU2D^2, na.rm = TRUE)
  energy_ok <- is.finite(row_energy) & row_energy > 0
  
  valid <- which(finite_ok & scatter_ok & energy_ok)
  if (!length(valid)) stop("No valid pixels after filtering.")
  
  X <- IFU2D[valid, , drop = FALSE]
  
  # ---- safe row scaling (as you wrote)
  safe_scale <- function(v) {
    ok <- is.finite(v)
    vv <- v[ok]
    if (length(vv) < 2) return(rep(0, length(v)))
    med <- stats::median(vv)
    sc  <- stats::mad(vv, center = med, constant = 1)
    if (!is.finite(sc) || sc == 0) sc <- stats::sd(vv)
    if (!is.finite(sc) || sc == 0) return(rep(0, length(v)))
    out <- (v - med) / sc
    out[!is.finite(out)] <- 0
    out
  }
  Xs <- t(apply(X, 1, safe_scale))
  
  # ---- PCA on scaled data (use irlba for big matrices)
  # installs: irlba
  pca <- irlba::prcomp_irlba(Xs, n = pca_k, center = FALSE, scale. = FALSE)
  Z <- pca$x  # (n_valid x pca_k)
  
  # ---- sample for Ward
  set.seed(42)
  m <- min(sample_n, nrow(Z))
  samp_idx <- sample.int(nrow(Z), m)
  Zs <- Z[samp_idx, , drop = FALSE]
  
  d <- dist(Zs)  # m*(m-1)/2 distances only
  hc <- fastcluster::hclust(d, method = "ward.D2")
  cl_s <- cutree(hc, k = Ncomp)
  
  # ---- centroids in PCA space (from sample)
  centers <- vapply(1:Ncomp, function(k) {
    colMeans(Zs[cl_s == k, , drop = FALSE])
  }, numeric(pca_k))
  centers <- t(centers)  # (Ncomp x pca_k)
  
  # ---- assign all valid pixels to nearest centroid (blockwise)
  # distance^2 = ||z||^2 + ||c||^2 - 2 z·c
  c2 <- rowSums(centers^2)
  z2 <- rowSums(Z^2)
  
  block <- 5000L
  assign <- integer(nrow(Z))
  for (i in seq(1, nrow(Z), by = block)) {
    j <- min(i + block - 1L, nrow(Z))
    ZZ <- Z[i:j, , drop = FALSE]
    # crossprod gives (block x Ncomp)
    dots <- ZZ %*% t(centers)
    # broadcast: z2 + c2 - 2*dots
    D2 <- (z2[i:j] + rep(c2, each = j - i + 1L)) - 2 * as.vector(dots)
    D2 <- matrix(D2, nrow = j - i + 1L, byrow = FALSE)
    assign[i:j] <- max.col(-D2)  # index of min distance
  }
  
  # ---- build map
  cluster_map <- matrix(NA_integer_, nrow = n_row, ncol = n_col)
  cluster_map[valid] <- assign
  
  list(cluster_map = cluster_map, valid = valid, pca = pca, centers = centers)
}
