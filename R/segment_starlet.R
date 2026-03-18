#' Segment a Cube After Applying a Sagui-style Starlet Mask
#'
#' This adds the photometric masking layer used by \code{sagui} before running
#' Capivara's spectral clustering. The cube is collapsed to white light,
#' decomposed with the starlet transform, reconstructed from selected scales,
#' masked spatially, and then segmented with \code{\link{segment_masked}}.
#'
#' @param input A FITS-like object with an \code{imDat} cube, or a raw
#'   3-D numeric array.
#' @param Ncomp Integer, the number of clusters to form.
#' @param redshift Numeric redshift placeholder kept for API compatibility.
#' @param scale_fn Row-wise scaling function used during segmentation.
#' @param collapse_fn Function used to build the white-light image.
#' @param starlet_J Number of starlet scales.
#' @param starlet_scales Integer vector of scales kept in the reconstruction.
#' @param include_coarse Logical; include the coarse starlet plane.
#' @param denoise_k Optional denoising threshold in MAD units.
#' @param mode Thresholding mode for starlet reconstruction.
#' @param positive_only Logical; keep only positive reconstructed values.
#' @param mask_mode Either \code{"na"} or \code{"zero"} for masked spaxels.
#'
#' @return A standard segmentation result plus the starlet products used to
#'   build the mask.
#' @export
segment_starlet <- function(input,
                            Ncomp = 5,
                            redshift = 0,
                            scale_fn = median_scale,
                            collapse_fn = collapse_white_light,
                            starlet_J = 5,
                            starlet_scales = 2:5,
                            include_coarse = FALSE,
                            denoise_k = 0,
                            mode = c("soft", "hard"),
                            positive_only = TRUE,
                            mask_mode = c("na", "zero")) {
  mode <- match.arg(mode)
  mask_mode <- match.arg(mask_mode)

  cubedat <- .as_cubedat(input)
  mask_info <- build_starlet_mask(
    input = cubedat,
    collapse_fn = collapse_fn,
    starlet_J = starlet_J,
    starlet_scales = starlet_scales,
    include_coarse = include_coarse,
    denoise_k = denoise_k,
    mode = mode,
    positive_only = positive_only
  )
  if (!any(mask_info$mask, na.rm = TRUE)) {
    stop("The starlet mask is empty. Try different `starlet_scales`, `denoise_k`, or `positive_only = FALSE`.")
  }

  masked_cube <- mask_cube(cubedat$imDat, mask_info$mask, mode = mask_mode)
  masked_input <- cubedat
  masked_input$imDat <- masked_cube

  out <- segment_masked(
    input = masked_input,
    Ncomp = Ncomp,
    redshift = redshift,
    scale_fn = scale_fn
  )

  c(out, list(
    mask = mask_info$mask,
    collapsed = mask_info$collapsed,
    starlet = list(
      decomposition = mask_info$decomposition,
      reconstruction = mask_info$reconstruction
    ),
    masked_cube = masked_cube
  ))
}
