#' Run a batch of Capivara-native barred MaNGA kinematic analyses
#'
#' @param batch_table_file CSV table with one row per galaxy.
#' @param output_root Root output directory for the master summary.
#' @return The master summary data frame.
#' @export
run_batch_barred_manga <- function(batch_table_file,
                                   output_root = "outputs") {
  .capivara_require(c("readr", "purrr"))
  batch <- readr::read_csv(batch_table_file, show_col_types = FALSE)
  .capivara_dir_create(output_root)

  rows <- vector("list", nrow(batch))
  for (i in seq_len(nrow(batch))) {
    row <- as.list(batch[i, , drop = FALSE])
    plateifu <- as.character(row$plateifu %||% paste0("row_", i))
    row$output_dir <- row$output_dir %||% file.path(output_root, plateifu)
    row$geometry <- list(
      x0 = row$x0 %||% NULL,
      y0 = row$y0 %||% NULL,
      pa_deg = row$pa_deg %||% NULL,
      inc_deg = row$inc_deg %||% NULL,
      vsys = row$vsys %||% NULL
    )
    row$bar <- list(
      bar_pa_deg = row$bar_pa_deg %||% NULL,
      bar_radius_pix = row$bar_radius_pix %||% NULL
    )
    out <- try(run_one_galaxy(row), silent = TRUE)
    if (inherits(out, "try-error")) {
      rows[[i]] <- data.frame(
        plateifu = plateifu,
        fit_status = "failed",
        error_message = conditionMessage(attr(out, "condition")),
        stringsAsFactors = FALSE
      )
    } else {
      rows[[i]] <- out$diagnostics$summary
      rows[[i]]$error_message <- NA_character_
    }
  }

  master <- dplyr::bind_rows(rows)
  master_file <- file.path(output_root, "barred_manga_capivara_kinematic_summary.csv")
  readr::write_csv(master, master_file)
  master
}

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L || all(is.na(x))) y else x
}
