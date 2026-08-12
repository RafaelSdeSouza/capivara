#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
site_dir <- if (length(args)) args[[1]] else "docs"

redirects <- c(
  "about.html" = "index.html#about",
  "news.html" = "news/index.html"
)

escape_html <- function(x) {
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  gsub('"', "&quot;", x, fixed = TRUE)
}

write_redirect <- function(from, to) {
  output <- file.path(site_dir, from)
  dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
  target <- escape_html(to)
  js_target <- encodeString(to, quote = '"')

  html <- c(
    "<!doctype html>",
    '<html lang="en">',
    "<head>",
    '<meta charset="utf-8">',
    '<meta name="viewport" content="width=device-width, initial-scale=1">',
    sprintf('<meta http-equiv="refresh" content="0; url=%s">', target),
    sprintf("<title>Moved to %s</title>", target),
    sprintf(
      "<script>var target=%s;var parts=target.split('#');var hash=window.location.hash||((parts.length>1)?'#'+parts.slice(1).join('#'):'');window.location.replace(parts[0]+window.location.search+hash);</script>",
      js_target
    ),
    "</head>",
    "<body>",
    "<main>",
    "<h1>This page has moved</h1>",
    sprintf('<p><a href="%s">Continue to the current documentation</a>.</p>', target),
    "</main>",
    "</body>",
    "</html>"
  )

  writeLines(html, output, useBytes = TRUE)
  message(from, " -> ", to)
}

if (!dir.exists(site_dir)) {
  stop("Site directory does not exist: ", site_dir)
}

Map(write_redirect, names(redirects), unname(redirects))
