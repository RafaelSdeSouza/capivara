#' Cluster a 2D Representation of a Data Cube
#'
#' This function processes a data cube (such as an IFU cube) by flattening it into
#' rows (spatial pixels) and columns (spectral variables), scaling each row, computing
#' pairwise distances using Capivara's internal distance helper, and performing hierarchical
#' clustering. The resulting clusters are rearranged back into a 2D grid consistent
#' with the original spatial dimensions. The function also retains the original data cube
#' for reference and post-processing.
#'
#' Missing spectral channels are handled automatically: non-finite values produced
#' during row-wise scaling are replaced with zero before distances are computed.
#' This keeps the standard exact workflow usable for masked cubes without a
#' separate public entry point.
#'
#' @param input A FITS object representing the input data cube. Typically, this is an IFU data cube.
#' @param Ncomp Integer, the number of clusters to form. Defaults to `15`.
#' @param redshift Systemic redshift. Required explicitly for rest-frame selection.
#' @param scale_fn A function used to scale each row of the 2D representation of the data cube.
#'   Defaults to \code{\link[base]{scale}}. If you have a custom scaling function, pass it here.
#' @param target_snr Optional minimum accepted SNR per cluster. When supplied,
#'   Capivara chooses the largest number of clusters whose minimum cluster SNR
#'   remains above this threshold.
#' @param var_cube Optional variance cube matching the input cube. Used only
#'   when \code{target_snr} is supplied.
#' @param k_values Optional candidate cluster counts tested when
#'   \code{target_snr} is supplied.
#' @param wavelength_range Optional wavelength interval used to compute SNR when
#'   \code{target_snr} is supplied. This SNR interval is always in native
#'   input coordinates, independently of `feature_wavelength_frame`.
#' @param wavelength_frame Frame of the input wavelength axis: `observed` or
#'   `rest`. May instead be supplied as `input$wavelength_frame`.
#' @param feature_wavelength_frame Frame of the requested feature interval:
#'   `observed` or `rest`. Required when bounds are supplied. Rest selection
#'   requires an explicit valid redshift, including for a rest-frame input.
#' @param feature_wavelength_range Optional wavelength interval used to select
#'   the spectral channels used for clustering. Selection includes native
#'   channels inside the closed interval; no resampling is performed.
#'   `wavelength_provenance` records frames, redshift, actual bounds and indices.
#'   \code{original_cube} remains the full input cube so downstream summed
#'   spectra are still flux-preserving across the full spectral axis.
#' @param snr_stat Either integrated SNR or median per-wavelength SNR when
#'   \code{target_snr} is supplied.
#' @param variance_inflation Multiplicative factor applied to propagated
#'   variances when \code{target_snr} is supplied.
#' @param use_starlet_mask Logical; if \code{TRUE}, build a Sagui-style
#'   support mask before clustering.
#' @param support_method Foreground support builder used when
#'   \code{use_starlet_mask = TRUE}. Options are \code{"starlet"} and
#'   \code{"adaptive"}.
#' @param support_args Optional named list passed to the selected support
#'   builder. For \code{"adaptive"}, arguments are passed to
#'   \code{\link{build_adaptive_support}}.
#' @param collapse_fn Function used to collapse the cube to white light when
#'   \code{use_starlet_mask = TRUE} and \code{support_method = "starlet"}.
#' @param starlet_J Number of starlet scales when \code{use_starlet_mask = TRUE}.
#' @param starlet_scales Integer vector of scales kept in the reconstruction
#'   when \code{use_starlet_mask = TRUE}.
#' @param include_coarse Logical; include the coarse starlet plane when
#'   \code{use_starlet_mask = TRUE}.
#' @param denoise_k Optional denoising threshold in MAD units when
#'   \code{use_starlet_mask = TRUE}.
#' @param starlet_mode Thresholding mode for the starlet reconstruction when
#'   \code{use_starlet_mask = TRUE}.
#' @param positive_only Logical; keep only positive reconstructed values when
#'   \code{use_starlet_mask = TRUE}.
#' @param mask_mode Either \code{"na"} or \code{"zero"} for masked spaxels
#'   when \code{use_starlet_mask = TRUE}.
#' @param support Optional `capivara_support` object. With an explicit
#'   representation, omission selects that representation's declared validity
#'   and eligibility domain; a supplied object is an additional explicit
#'   spatial restriction. Historical `representation = NULL` calls retain their
#'   original support behavior. Support cannot be combined with legacy starlet
#'   masking.
#' @param representation Optional explicit `capivara_representation`. When
#'   supplied, CAPIVARA uses the observed-entry spatial hierarchy; the default
#'   `NULL` retains the historical median-centred all-pairs implementation.
#' @param sample_validity Optional logical row-column-channel validity array.
#'   A support object's quality record takes precedence when this is omitted.
#'
#' @details
#' Steps performed by the function:
#' \enumerate{
#'   \item Reads the input FITS data cube.
#'   \item Converts the cube into a 2D matrix (spatial pixels x spectral variables).
#'   \item Scales the data row-wise using \code{scale_fn}.
#'   \item Computes pairwise distances between rows using an internal distance helper.
#'   \item Performs hierarchical clustering using Ward's D2 method via \code{\link[fastcluster]{hclust}}.
#'   \item Cuts the dendrogram into \code{Ncomp} clusters and reshapes the results into a 2D cluster map,
#'   or, when \code{target_snr} is supplied, chooses the largest cut whose minimum cluster SNR
#'   remains above the requested threshold.
#'   \item Calculates the signal-to-noise ratio (SNR) for each cluster.
#' }
#'
#' @return A list containing:
#' \item{cluster_map}{A \code{n_rows x n_cols} matrix of cluster assignments (integers), where each element corresponds to a spatial pixel in the original cube layout.}
#' \item{header}{The header metadata from the input FITS file.}
#' \item{axDat}{Axis information (e.g., spatial and spectral axes) from the input FITS file.}
#' \item{cluster_snr}{A numeric vector containing the signal-to-noise ratio (SNR) for each cluster.}
#' \item{original_cube}{The original FITS data cube as input to the function, for reference and post-processing.}
#' \item{starlet_info}{When \code{use_starlet_mask = TRUE}, a list containing
#' the spatial mask, white-light image, starlet decomposition,
#' reconstruction, and masked cube used before clustering.}
#'
#' This process is often used in IFU data analysis, where clustering is applied to grouped
#' spectral profiles of spatial pixels to identify regions with similar characteristics.
#'
#' @seealso \code{\link{segment_large}}, \code{\link{build_starlet_mask}},
#'   \code{\link[fastcluster]{hclust}},
#'   \code{\link[stats]{cutree}}
#'
#' @examples
#' input_cube <- list(imDat = array(runif(5 * 5 * 12), dim = c(5, 5, 12)))
#' clustering_result <- segment(input = input_cube, Ncomp = 5)
#' cluster_map <- clustering_result$cluster_map
#' original_cube <- clustering_result$original_cube
#'
#' @export
segment <- function(input,
                    Ncomp = 15,
                    redshift = NA_real_,
                    scale_fn = median_scale,
                    target_snr = NULL,
                    var_cube = NULL,
                    k_values = NULL,
                    wavelength_range = NULL,
                    feature_wavelength_range = NULL,
                    snr_stat = c("integrated", "median_per_wavelength"),
                    variance_inflation = 1,
                    use_starlet_mask = FALSE,
                    support_method = c("starlet", "adaptive"),
                    support_args = list(),
                    collapse_fn = collapse_white_light,
                    starlet_J = 5,
                    starlet_scales = 2:5,
                    include_coarse = FALSE,
                    denoise_k = 0,
                    starlet_mode = c("soft", "hard"),
                    positive_only = TRUE,
                    mask_mode = c("na", "zero"),
                    wavelength_frame = NULL,
                    feature_wavelength_frame = NULL,
                    support = NULL,
                    representation = NULL,
                    sample_validity = NULL) {
  scale_fn_supplied <- !missing(scale_fn)
  starlet_mode <- match.arg(starlet_mode)
  mask_mode <- match.arg(mask_mode)
  support_method <- match.arg(support_method)
  if (!is.null(representation)) {
    if (!is.null(target_snr)) stop("Representation-specific hierarchy does not yet implement target-S/N cutting; use a fixed chronological `Ncomp`.", call. = FALSE)
    if (scale_fn_supplied) stop("An explicit representation owns its transformation; do not also supply `scale_fn`.", call. = FALSE)
    if (!is.null(feature_wavelength_range) || !is.null(feature_wavelength_frame)) stop("The representation declares its coordinate domain; do not also supply legacy feature bounds.", call. = FALSE)
    if (isTRUE(use_starlet_mask)) stop("Build and supply one explicit support object before representation-specific segmentation.", call. = FALSE)
    return(.segment_semantic(input, Ncomp, redshift, var_cube, representation,
                             sample_validity, support, wavelength_frame,
                             return_details = TRUE))
  }
  if (!is.null(support) && isTRUE(use_starlet_mask)) {
    stop("Supply a `capivara_support` object or request legacy starlet masking, not both.")
  }

  raw_input <- .as_cubedat(input)

  starlet_prep <- .apply_starlet_support(
    input = input,
    use_starlet_mask = use_starlet_mask,
    support_method = support_method,
    support_args = support_args,
    collapse_fn = collapse_fn,
    starlet_J = starlet_J,
    starlet_scales = starlet_scales,
    include_coarse = include_coarse,
    denoise_k = denoise_k,
    mode = starlet_mode,
    positive_only = positive_only,
    mask_mode = mask_mode
  )
  supported_input <- .apply_capivara_support(starlet_prep$input, support)
  full_input <- if (is.null(support)) supported_input else raw_input
  feature_subset <- .subset_cubedat_wavelength_range(
    supported_input,
    feature_wavelength_range = feature_wavelength_range,
    wavelength_frame = wavelength_frame,
    feature_wavelength_frame = feature_wavelength_frame,
    redshift = redshift
  )
  input <- feature_subset$cubedat

  if (!is.null(var_cube) && !is.null(feature_wavelength_range)) {
    var_cube <- .subset_variance_channels(var_cube, full_input, feature_subset$wave_idx)
  }

  if (!is.null(target_snr)) {
    if (!missing(Ncomp)) {
      stop("Specify either `Ncomp` or `target_snr`, not both.")
    }

    out <- choose_ncomp_by_snr(
      input = input,
      target_snr = target_snr,
      var_cube = var_cube,
      k_values = k_values,
      wavelength_range = wavelength_range,
      redshift = redshift,
      scale_fn = scale_fn,
      snr_stat = snr_stat,
      variance_inflation = variance_inflation,
      na_safe = TRUE
    )
  } else {
    out <- .segment_core(
      input = input,
      Ncomp = Ncomp,
      redshift = redshift,
      scale_fn = scale_fn,
      na_to_zero = TRUE
    )
  }

  if (!is.null(starlet_prep$starlet_info)) {
    out$starlet_info <- starlet_prep$starlet_info
  }
  if (!is.null(starlet_prep$support_info)) {
    out$support_info <- starlet_prep$support_info
  }
  if (!is.null(support)) {
    warning("`support` without `representation` restricts spatial enrollment but retains historical spectral scaling and missing-value behavior.", call. = FALSE)
    out$support <- support
    out$support_id <- support$support_id
    out$support_provenance <- support$provenance
    out$cluster_map[!support$analysis_mask] <- NA_integer_
  }

  out$original_cube <- full_input
  out$header <- full_input$hdr
  out$axDat <- full_input$axDat

  out$wavelength_provenance <- feature_subset$provenance

  if (!is.null(feature_wavelength_range)) {
    out$feature_wavelength_range <- feature_wavelength_range
    out$feature_wavelength_index <- feature_subset$wave_idx
    out$feature_wavelengths <- feature_subset$selected_wavelengths
  }

  out <- .capivara_attach_legacy_provenance(
    out, support, starlet_prep$starlet_info, starlet_prep$support_info,
    validity_rule = "positive finite integrated signal on the historical input path",
    missing_value_rule = "historical row scaling with non-finite entries replaced by zero"
  )
  .capivara_attach_spatial_products(out, full_input)
}
