#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
site_dir <- normalizePath(if (length(args)) args[[1]] else "docs", mustWork = TRUE)

if (!requireNamespace("xml2", quietly = TRUE)) {
  stop("The xml2 package is required.")
}

html_files <- list.files(site_dir, pattern = "[.]html$", recursive = TRUE, full.names = TRUE)
if (!length(html_files)) stop("No HTML files found in ", site_dir)

normalise_target <- function(source, href) {
  href <- sub("[?].*$", "", href)
  path <- sub("#.*$", "", href)
  fragment <- if (grepl("#", href, fixed = TRUE)) sub("^[^#]*#", "", href) else ""

  if (!nzchar(path)) {
    target <- source
  } else if (startsWith(path, "/")) {
    target <- file.path(site_dir, sub("^/+", "", path))
  } else {
    target <- file.path(dirname(source), utils::URLdecode(path))
  }

  if (dir.exists(target)) target <- file.path(target, "index.html")
  if (!file.exists(target) && !grepl("[.][A-Za-z0-9]+$", target)) {
    html_target <- paste0(target, ".html")
    if (file.exists(html_target)) target <- html_target
  }

  list(path = normalizePath(target, mustWork = FALSE), fragment = utils::URLdecode(fragment))
}

errors <- character()
for (source in html_files) {
  doc <- tryCatch(xml2::read_html(source), error = function(e) e)
  if (inherits(doc, "error")) {
    errors <- c(errors, sprintf("Unreadable HTML: %s (%s)", source, conditionMessage(doc)))
    next
  }

  is_redirect <- length(xml2::xml_find_all(
    doc,
    ".//meta[translate(@http-equiv, 'ABCDEFGHIJKLMNOPQRSTUVWXYZ', 'abcdefghijklmnopqrstuvwxyz')='refresh']"
  )) > 0L
  h1_count <- length(xml2::xml_find_all(doc, ".//h1"))
  if (!is_redirect && h1_count != 1L) {
    errors <- c(errors, sprintf("Expected one H1, found %d: %s", h1_count, source))
  }

  images <- xml2::xml_find_all(doc, ".//main//img")
  if (length(images)) {
    alt <- trimws(xml2::xml_attr(images, "alt"))
    missing_alt <- is.na(alt) | !nzchar(alt)
    if (any(missing_alt)) {
      errors <- c(errors, sprintf("Missing main-content image alt text: %s", source))
    }
  }

  links <- unique(xml2::xml_attr(xml2::xml_find_all(doc, ".//a[@href]"), "href"))
  links <- links[!is.na(links) & nzchar(links)]
  links <- links[!grepl("^(https?:|//|mailto:|tel:|javascript:|data:)", links, ignore.case = TRUE)]

  for (href in links) {
    target <- normalise_target(source, href)
    if (!file.exists(target$path)) {
      errors <- c(errors, sprintf("Missing target: %s -> %s", source, href))
      next
    }
    if (nzchar(target$fragment) && grepl("[.]html$", target$path, ignore.case = TRUE)) {
      target_doc <- tryCatch(xml2::read_html(target$path), error = function(e) NULL)
      if (is.null(target_doc)) next
      ids <- c(
        xml2::xml_attr(xml2::xml_find_all(target_doc, ".//*[@id]"), "id"),
        xml2::xml_attr(xml2::xml_find_all(target_doc, ".//a[@name]"), "name")
      )
      if (!target$fragment %in% ids) {
        errors <- c(errors, sprintf("Missing fragment: %s -> %s", source, href))
      }
    }
  }
}

if (length(errors)) {
  cat(paste(unique(errors), collapse = "\n"), "\n")
  quit(status = 1L)
}

cat(sprintf("PASS: checked %d HTML files for internal targets, fragments, H1s, and main-content alt text.\n", length(html_files)))
