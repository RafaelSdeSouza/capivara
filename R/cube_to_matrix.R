#' Transform a 3D IFU Cube into a 2D Matrix
#'
#' This function reshapes a 3D integral field unit (IFU) data cube into a 2D matrix.
#' The input cube is expected to have dimensions corresponding to two spatial axes
#' and one spectral axis (e.g., \code{[x, y, wavelength]}). The output is a matrix
#' where each row represents a single spatial pixel, and each column corresponds to
#' flux values across wavelengths.
#'
#' @param x An object containing IFU data. Must include an element \code{imDat},
#'   which is a 3D numeric array of dimensions \eqn{(n\_row, n\_col, n\_wave)}.
#'
#' @details
#' This function takes the input cube \code{x$imDat} and flattens the spatial dimensions
#' (x and y) into one dimension. The resulting matrix has \eqn{n\_row * n\_col} rows,
#' each corresponding to a single spatial pixel (with the pixel indexing following the
#' array's column-major order used by R). The \eqn{n\_wave} dimension is left intact as
#' the number of columns, preserving the spectral information. If the spatial or spectral
#' arrangement in your data differ, you may need to permute or rotate the data before
#' calling this function.
#'
#' @return A 2D numeric matrix of dimensions \code{(n\_row * n\_col) x n\_wave}.
#'
#' @seealso
#' Related functions in IFU processing packages that handle cube I/O, slicing, and
#' spectral extraction.
#'
#' @examples
#' # Create a small synthetic IFU cube: 10x10 pixels, 50 wavelengths
#' sim_cube <- array(runif(10 * 10 * 50), dim = c(10, 10, 50))
#' x <- list(imDat = sim_cube)
#'
#' # Convert the 3D cube into a 2D matrix
#' result <- cube_to_matrix(x)
#' dim(result)  # Should be c(100, 50)
#'
#' @export
cube_to_matrix <- function(x) {
  IFS_gal <- as.array(x$imDat)

  # Check that IFS_gal is 3D
  if (length(dim(IFS_gal)) != 3) {
    stop("Error: 'x$imDat' must be a 3D array with dimensions (n_row, n_col, n_wave).")
  }

  n_row <- dim(IFS_gal)[1]
  n_col <- dim(IFS_gal)[2]
  n_wave <- dim(IFS_gal)[3]

  # Reshape the array into a matrix: rows = pixels, columns = wavelengths
  dim(IFS_gal) <- c(n_row * n_col, n_wave)

  return(IFS_gal)
}

