#' Convert optical laboratory wavelength between air and vacuum
#'
#' Uses the standard-air refractive index of Morton (1991, ApJS 77, 119), as tabulated by SDSS DR5,
#' evaluated at vacuum wavelength. The inverse is solved iteratively; it is
#' not the approximation obtained by evaluating the index at air wavelength.
#' This changes coordinates only, not sampled flux densities or LSF widths.
#' Apply it to laboratory wavelengths before applying a physical redshift.
#' @param wavelength Numeric wavelengths in Angstrom, between 2000 and 20000.
#' @param from,to Explicit input and output medium, `air` or `vacuum`.
#' @return Numeric converted wavelengths in Angstrom.
#' @export
convert_wavelength_medium <- function(wavelength, from, to) {
  from <- match.arg(from, c("air", "vacuum"))
  to <- match.arg(to, c("air", "vacuum"))
  if (!is.numeric(wavelength) || any(!is.finite(wavelength)) ||
      any(wavelength < 2000 | wavelength > 20000)) {
    stop("Air/vacuum conversion requires finite optical wavelengths (2000--20000 Angstrom).")
  }
  if (from == to) return(wavelength)
  index <- function(vac) {
    1 + 2.735182e-4 + 131.4182/vac^2 + 2.76249e8/vac^4
  }
  if (from == "vacuum") return(wavelength/index(wavelength))
  vac <- wavelength
  for (i in seq_len(8)) vac <- wavelength * index(vac)
  vac
}

.capivara_emission_line_table <- function(wavelength_medium = c("vacuum", "air")) {
  medium <- match.arg(wavelength_medium)
  tab <- data.frame(
    name = c("oii3726", "oii3729", "neiii3869", "hdelta", "hgamma", "hbeta",
      "oiii4959", "oiii5007", "oi6300", "halpha", "nii6548", "nii6583", "sii6716", "sii6731"),
    label = c("[O II] 3726", "[O II] 3729", "[Ne III] 3869", "Hdelta", "Hgamma", "Hbeta",
      "[O III] 4959", "[O III] 5007", "[O I] 6300", "Halpha", "[N II] 6548", "[N II] 6583", "[S II] 6716", "[S II] 6731"),
    rest_wavelength = c(3727.0920, 3729.8750, 3869.8600, 4102.8922, 4341.6837, 4862.6830,
      4960.2950, 5008.2400, 6302.0460, 6564.6080, 6549.8600, 6585.2700, 6718.2950, 6732.6740),
    family = c("blue", "blue", "blue", "balmer", "balmer", "balmer", "agn", "agn", "agn", "balmer", rep("agn", 4)),
    wavelength_medium = medium,
    reference = "https://www.sdss4.org/dr17/manga/manga-analysis-pipeline/ (Table 2)",
    registry_version = "SDSS-DR17-DAP-Table2",
    stringsAsFactors = FALSE)
  tab$source_vacuum_wavelength <- tab$rest_wavelength
  tab$conversion <- "none"
  if (medium == "air") {
    tab$rest_wavelength <- convert_wavelength_medium(tab$rest_wavelength, "vacuum", "air")
    tab$conversion <- "SDSS DR5 / Morton 1991; https://classic.sdss.org/dr5/products/spectra/vacwavelength.php"
  }
  tab
}
