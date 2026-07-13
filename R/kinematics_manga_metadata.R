#' Infer a MaNGA plate-IFU identifier from text
#'
#' @param x File name, path, header value, or plate-IFU-like string.
#' @return A character scalar like `"8078-12703"`, or `NA_character_`.
#' @noRd
infer_manga_plateifu <- function(x) {
  if (is.null(x) || !length(x) || is.na(x[[1]])) {
    return(NA_character_)
  }
  x <- basename(trimws(as.character(x[[1]])))
  hit <- regmatches(x, regexpr("[0-9]{4,5}-[0-9]{3,5}", x))
  if (length(hit) && nzchar(hit[[1]])) hit[[1]] else NA_character_
}

.manga_metadata_file <- function(metadata_file = NULL) {
  file_name <- "manga_drpall_v3_1_1_redshifts.csv.gz"
  candidates <- character()
  if (!is.null(metadata_file) && length(metadata_file) && !is.na(metadata_file[[1]]) && nzchar(metadata_file[[1]])) {
    candidates <- c(candidates, metadata_file[[1]])
  }
  env_file <- Sys.getenv("CAPIVARA_MANGA_METADATA", unset = "")
  if (nzchar(env_file)) {
    candidates <- c(candidates, env_file)
  }
  installed <- system.file("extdata", file_name, package = "capivara", mustWork = FALSE)
  if (nzchar(installed)) {
    candidates <- c(candidates, installed)
  }
  repo <- Sys.getenv("CAPIVARA_REPO_ROOT", unset = "")
  if (nzchar(repo)) {
    candidates <- c(candidates, file.path(repo, "inst", "extdata", file_name))
  }
  candidates <- c(
    candidates,
    file.path(getwd(), "inst", "extdata", file_name)
  )
  hit <- candidates[file.exists(candidates)][1]
  if (length(hit) && !is.na(hit)) {
    normalizePath(hit, mustWork = TRUE)
  } else {
    stop("Could not find MaNGA redshift metadata table. Set CAPIVARA_MANGA_METADATA to a DRPall-derived CSV.", call. = FALSE)
  }
}

.manga_clean_redshift <- function(x) {
  z <- suppressWarnings(as.numeric(x))
  z[!is.finite(z) | z <= -0.01 | z >= 2] <- NA_real_
  z
}

#' Read the bundled MaNGA DRPall redshift table
#'
#' @param metadata_file Optional CSV/CSV.GZ file with `plateifu` and `redshift`.
#' @return A data frame with MaNGA identifiers and redshifts.
#' @noRd
read_manga_redshift_table <- function(metadata_file = NULL) {
  path <- .manga_metadata_file(metadata_file)
  tab <- utils::read.csv(path, stringsAsFactors = FALSE)
  tab$plateifu <- trimws(as.character(tab$plateifu))
  tab$mangaid <- trimws(as.character(tab$mangaid))
  tab$redshift <- .manga_clean_redshift(tab$redshift)
  if ("nsa_z" %in% names(tab)) {
    tab$nsa_z <- .manga_clean_redshift(tab$nsa_z)
  }
  tab
}

#' Look up local MaNGA metadata
#'
#' @param plateifu Plate-IFU identifier, e.g. `"8078-12703"`.
#' @param mangaid Optional MaNGA-ID fallback.
#' @param metadata_file Optional local metadata CSV/CSV.GZ.
#' @return A one-row data frame.
#' @noRd
manga_metadata <- function(plateifu = NULL, mangaid = NULL, metadata_file = NULL) {
  tab <- read_manga_redshift_table(metadata_file)
  keep <- rep(FALSE, nrow(tab))
  if (!is.null(plateifu) && length(plateifu) && !is.na(plateifu[[1]]) && nzchar(plateifu[[1]])) {
    target <- tolower(infer_manga_plateifu(plateifu[[1]]))
    keep <- keep | tolower(tab$plateifu) == target
  }
  if (!is.null(mangaid) && length(mangaid) && !is.na(mangaid[[1]]) && nzchar(mangaid[[1]])) {
    keep <- keep | tolower(tab$mangaid) == tolower(trimws(as.character(mangaid[[1]])))
  }
  out <- tab[keep, , drop = FALSE]
  if (!nrow(out)) {
    stop("No MaNGA metadata row found for plateifu/mangaid.", call. = FALSE)
  }
  out[1, , drop = FALSE]
}

#' Look up a local MaNGA redshift
#'
#' @param plateifu Plate-IFU identifier or cube filename.
#' @param mangaid Optional MaNGA-ID fallback.
#' @param metadata_file Optional local metadata CSV/CSV.GZ.
#' @return Numeric redshift.
#' @noRd
manga_redshift <- function(plateifu = NULL, mangaid = NULL, metadata_file = NULL) {
  row <- manga_metadata(plateifu = plateifu, mangaid = mangaid, metadata_file = metadata_file)
  z <- row$redshift[[1]]
  if (!is.finite(z) && "nsa_z" %in% names(row)) {
    z <- row$nsa_z[[1]]
  }
  if (!is.finite(z)) {
    stop("Metadata row exists, but no finite redshift is available.", call. = FALSE)
  }
  z
}

.manga_header_value <- function(hdr, key, default = NA_character_) {
  idx <- which(hdr == key)[1]
  if (!length(idx) || is.na(idx) || idx >= length(hdr)) {
    return(default)
  }
  value <- hdr[idx + 1L]
  if (is.na(value) || !nzchar(value)) default else value
}

#' Read basic MaNGA cube metadata from a LOGCUBE header
#'
#' @param cube_path Path to a MaNGA LOGCUBE.
#' @return A list with `plateifu`, `mangaid`, and `redshift`.
#' @noRd
read_manga_cube_metadata <- function(cube_path) {
  plateifu <- infer_manga_plateifu(cube_path)
  mangaid <- NA_character_
  redshift <- NA_real_
  if (file.exists(cube_path) && requireNamespace("FITSio", quietly = TRUE)) {
    hdr <- tryCatch(FITSio::readFITS(cube_path, hdu = 1L, maxLines = 1)$hdr, error = function(e) NULL)
    if (!is.null(hdr)) {
      plateifu <- infer_manga_plateifu(.manga_header_value(hdr, "PLATEIFU", plateifu))
      mangaid <- trimws(.manga_header_value(hdr, "MANGAID", NA_character_))
      for (key in c("NSA_Z", "REDSHIFT", "Z", "Z_SPEC", "SPECZ", "DRP_Z", "OBJZ")) {
        z <- .manga_clean_redshift(.manga_header_value(hdr, key, NA_character_))
        if (is.finite(z)) {
          redshift <- z
          break
        }
      }
    }
  }
  list(plateifu = plateifu, mangaid = mangaid, redshift = redshift)
}

#' Resolve the redshift for a MaNGA cube
#'
#' Resolution order is manual value, FITS header, then bundled/local DRPall
#' metadata.
#'
#' @param cube_path Path to a MaNGA LOGCUBE.
#' @param redshift Optional manual redshift.
#' @param metadata_file Optional local metadata CSV/CSV.GZ.
#' @return A list with `redshift`, `source`, `plateifu`, and `mangaid`.
#' @noRd
resolve_manga_redshift <- function(cube_path, redshift = NA_real_, metadata_file = NULL) {
  manual <- .manga_clean_redshift(redshift)
  if (is.finite(manual)) {
    meta <- read_manga_cube_metadata(cube_path)
    return(list(redshift = manual, source = "manual", plateifu = meta$plateifu, mangaid = meta$mangaid))
  }
  meta <- read_manga_cube_metadata(cube_path)
  if (is.finite(meta$redshift)) {
    return(list(redshift = meta$redshift, source = "fits_header", plateifu = meta$plateifu, mangaid = meta$mangaid))
  }
  z <- manga_redshift(plateifu = meta$plateifu, mangaid = meta$mangaid, metadata_file = metadata_file)
  list(redshift = z, source = "local_drpall_v3_1_1", plateifu = meta$plateifu, mangaid = meta$mangaid)
}
