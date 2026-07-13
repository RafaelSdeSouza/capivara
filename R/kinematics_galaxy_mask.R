#' Clean a logical IFU footprint into a galaxy mask
#'
#' This is intentionally conservative: it keeps the largest connected valid
#' footprint, fills enclosed holes, and can close small one-spaxel gaps. It is
#' meant to define where the velocity model is allowed to fit; starlet support
#' should usually be used as a soft weight on top of this mask, not as the
#' footprint itself.
#'
#' @param valid_mask Logical matrix of valid science spaxels.
#' @param close_iterations Number of dilate/erode closing passes.
#' @param fill_holes Fill enclosed holes after component selection.
#' @param preserve_input If TRUE, never remove pixels present in `valid_mask`.
#' @param connectivity Either 4 or 8.
#' @return Logical matrix with the same dimensions as `valid_mask`.
#' @export
better_galaxy_mask <- function(valid_mask,
                               close_iterations = 1L,
                               fill_holes = TRUE,
                               preserve_input = TRUE,
                               connectivity = 8L) {
  if (!is.matrix(valid_mask)) {
    stop("`valid_mask` must be a logical matrix.", call. = FALSE)
  }
  input_mask <- isTRUE(valid_mask) | (is.logical(valid_mask) & valid_mask)
  mask <- input_mask
  mask[is.na(mask)] <- FALSE
  input_mask[is.na(input_mask)] <- FALSE
  mask <- .capivara_largest_component(mask, connectivity = connectivity)
  if (isTRUE(fill_holes)) {
    mask <- .capivara_fill_holes(mask, connectivity = connectivity)
  }
  close_iterations <- as.integer(close_iterations)
  if (is.finite(close_iterations) && close_iterations > 0L) {
    for (i in seq_len(close_iterations)) {
      mask <- .capivara_erode_mask(.capivara_dilate_mask(mask, connectivity), connectivity)
      if (isTRUE(fill_holes)) {
        mask <- .capivara_fill_holes(mask, connectivity = connectivity)
      }
    }
  }
  if (isTRUE(preserve_input)) {
    mask <- mask | input_mask
  }
  mask
}

.capivara_neighbor_offsets <- function(connectivity = 8L) {
  connectivity <- as.integer(connectivity)
  if (connectivity == 4L) {
    return(rbind(c(-1L, 0L), c(1L, 0L), c(0L, -1L), c(0L, 1L)))
  }
  if (connectivity != 8L) {
    stop("`connectivity` must be 4 or 8.", call. = FALSE)
  }
  out <- as.matrix(expand.grid(dr = -1L:1L, dc = -1L:1L))
  out[rowSums(abs(out)) > 0, , drop = FALSE]
}

.capivara_shift_mask <- function(mask, dr, dc, fill = FALSE) {
  nr <- nrow(mask)
  nc <- ncol(mask)
  out <- matrix(fill, nr, nc)
  r_src <- seq_len(nr)
  c_src <- seq_len(nc)
  r_dst <- r_src + dr
  c_dst <- c_src + dc
  keep_r <- r_dst >= 1L & r_dst <= nr
  keep_c <- c_dst >= 1L & c_dst <= nc
  out[r_dst[keep_r], c_dst[keep_c]] <- mask[r_src[keep_r], c_src[keep_c]]
  out
}

.capivara_dilate_mask <- function(mask, connectivity = 8L) {
  out <- mask
  offsets <- .capivara_neighbor_offsets(connectivity)
  for (i in seq_len(nrow(offsets))) {
    out <- out | .capivara_shift_mask(mask, offsets[i, 1], offsets[i, 2], fill = FALSE)
  }
  out
}

.capivara_erode_mask <- function(mask, connectivity = 8L) {
  out <- mask
  offsets <- .capivara_neighbor_offsets(connectivity)
  for (i in seq_len(nrow(offsets))) {
    out <- out & .capivara_shift_mask(mask, offsets[i, 1], offsets[i, 2], fill = FALSE)
  }
  out
}

.capivara_largest_component <- function(mask, connectivity = 8L) {
  nr <- nrow(mask)
  nc <- ncol(mask)
  visited <- matrix(FALSE, nr, nc)
  labels <- matrix(0L, nr, nc)
  offsets <- .capivara_neighbor_offsets(connectivity)
  best_label <- 0L
  best_size <- 0L
  label <- 0L
  seeds <- which(mask, arr.ind = TRUE)
  if (!nrow(seeds)) {
    return(mask)
  }

  for (s in seq_len(nrow(seeds))) {
    r0 <- seeds[s, 1]
    c0 <- seeds[s, 2]
    if (visited[r0, c0]) {
      next
    }
    label <- label + 1L
    queue_r <- integer(nr * nc)
    queue_c <- integer(nr * nc)
    head <- 1L
    tail <- 1L
    queue_r[tail] <- r0
    queue_c[tail] <- c0
    visited[r0, c0] <- TRUE
    size <- 0L

    while (head <= tail) {
      r <- queue_r[head]
      c <- queue_c[head]
      head <- head + 1L
      labels[r, c] <- label
      size <- size + 1L
      for (i in seq_len(nrow(offsets))) {
        rr <- r + offsets[i, 1]
        cc <- c + offsets[i, 2]
        if (rr < 1L || rr > nr || cc < 1L || cc > nc) {
          next
        }
        if (!visited[rr, cc] && mask[rr, cc]) {
          tail <- tail + 1L
          queue_r[tail] <- rr
          queue_c[tail] <- cc
          visited[rr, cc] <- TRUE
        }
      }
    }
    if (size > best_size) {
      best_size <- size
      best_label <- label
    }
  }
  labels == best_label
}

.capivara_fill_holes <- function(mask, connectivity = 8L) {
  nr <- nrow(mask)
  nc <- ncol(mask)
  background <- !mask
  outside <- matrix(FALSE, nr, nc)
  offsets <- .capivara_neighbor_offsets(connectivity)
  border <- unique(rbind(
    cbind(1L, seq_len(nc)),
    cbind(nr, seq_len(nc)),
    cbind(seq_len(nr), 1L),
    cbind(seq_len(nr), nc)
  ))
  seeds <- border[background[border], , drop = FALSE]
  if (!nrow(seeds)) {
    return(mask)
  }

  queue_r <- integer(nr * nc)
  queue_c <- integer(nr * nc)
  head <- 1L
  tail <- 0L
  for (s in seq_len(nrow(seeds))) {
    r <- seeds[s, 1]
    c <- seeds[s, 2]
    if (!outside[r, c]) {
      tail <- tail + 1L
      queue_r[tail] <- r
      queue_c[tail] <- c
      outside[r, c] <- TRUE
    }
  }
  while (head <= tail) {
    r <- queue_r[head]
    c <- queue_c[head]
    head <- head + 1L
    for (i in seq_len(nrow(offsets))) {
      rr <- r + offsets[i, 1]
      cc <- c + offsets[i, 2]
      if (rr < 1L || rr > nr || cc < 1L || cc > nc) {
        next
      }
      if (background[rr, cc] && !outside[rr, cc]) {
        tail <- tail + 1L
        queue_r[tail] <- rr
        queue_c[tail] <- cc
        outside[rr, cc] <- TRUE
      }
    }
  }
  mask | (background & !outside)
}
