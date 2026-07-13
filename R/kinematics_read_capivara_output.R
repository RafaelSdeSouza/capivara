#' Read Capivara segmentation products for kinematic modelling
#'
#' @param segmentation_file FITS or RDS file containing a 2-D integer
#'   segmentation map.
#' @param segment_table_file CSV table mapping segment labels to classes and
#'   fitting flags.
#' @return A list with `segmentation_map` and `segment_table`.
#' @export
read_capivara_output <- function(segmentation_file, segment_table_file) {
  .capivara_require(c("FITSio", "readr"))

  if (!file.exists(segmentation_file)) {
    stop("Segmentation file does not exist: ", segmentation_file, call. = FALSE)
  }
  if (!file.exists(segment_table_file)) {
    stop("Segment table file does not exist: ", segment_table_file, call. = FALSE)
  }

  if (grepl("\\.rds$", segmentation_file, ignore.case = TRUE)) {
    obj <- readRDS(segmentation_file)
    segmentation_map <- if (is.list(obj) && !is.null(obj$segmentation_map)) {
      obj$segmentation_map
    } else {
      obj
    }
  } else {
    fits_path <- .fits_temp_unzip(segmentation_file)
    on.exit(if (!identical(fits_path, segmentation_file)) unlink(fits_path), add = TRUE)
    segmentation_map <- FITSio::readFITS(fits_path, hdu = 1, maxLines = 20000)$imDat
  }

  segmentation_map <- matrix(as.integer(round(segmentation_map)), nrow = nrow(segmentation_map))
  segment_table <- readr::read_csv(segment_table_file, show_col_types = FALSE)

  required <- c("label", "class")
  missing <- setdiff(required, names(segment_table))
  if (length(missing)) {
    stop("Segment table is missing required column(s): ",
         paste(missing, collapse = ", "), call. = FALSE)
  }

  segment_table$label <- as.integer(segment_table$label)
  segment_table$class <- as.character(segment_table$class)
  if (!"use_for_disc_fit" %in% names(segment_table)) {
    segment_table$use_for_disc_fit <- tolower(segment_table$class) == "disc"
  }
  if (!"use_for_bar_diagnostics" %in% names(segment_table)) {
    segment_table$use_for_bar_diagnostics <- tolower(segment_table$class) %in% c("bar", "ring")
  }
  segment_table$use_for_disc_fit <- as.logical(segment_table$use_for_disc_fit)
  segment_table$use_for_bar_diagnostics <- as.logical(segment_table$use_for_bar_diagnostics)

  labels_in_map <- sort(unique(as.integer(stats::na.omit(as.vector(segmentation_map)))))
  labels_in_map <- labels_in_map[labels_in_map > 0L]
  missing_labels <- setdiff(labels_in_map, segment_table$label)
  if (length(missing_labels)) {
    stop(
      "Segmentation labels missing from segment table: ",
      paste(missing_labels, collapse = ", "),
      call. = FALSE
    )
  }

  list(segmentation_map = segmentation_map, segment_table = as.data.frame(segment_table))
}
