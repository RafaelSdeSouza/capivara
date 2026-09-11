#!/usr/bin/env Rscript
# Compare the installed API/function bodies and bundled scripts with a source tree.
args <- commandArgs(trailingOnly=TRUE)
if (!length(args)) stop("Usage: Rscript tools/check_source_install.R SOURCE [LIBRARY]")
repo <- normalizePath(args[1],mustWork=TRUE)
if (length(args)>1) .libPaths(c(normalizePath(args[2],mustWork=TRUE),.libPaths()))
desc <- read.dcf(file.path(repo,"DESCRIPTION"))
package <- unname(desc[1,"Package"])
ns <- asNamespace(package)
installed <- packageDescription(package)
stopifnot(identical(installed$Version,unname(desc[1,"Version"])))
source_functions <- list()
for (f in list.files(file.path(repo,"R"),pattern="[.]R$",full.names=TRUE)) {
  for (expr in parse(f,keep.source=FALSE)) {
    if (is.call(expr) && identical(expr[[1]],as.name("<-")) &&
        is.name(expr[[2]]) && is.call(expr[[3]]) &&
        identical(expr[[3]][[1]],as.name("function"))) {
      source_functions[[as.character(expr[[2]])]] <- eval(expr[[3]],envir=ns)
    }
  }
}
canonical <- function(x)paste(deparse(x,width.cutoff=500L,control="all"),collapse="\n")
bad <- names(source_functions)[!vapply(names(source_functions),function(n){
  f <- get0(n,envir=ns,inherits=FALSE)
  is.function(f) && identical(canonical(formals(f)),canonical(formals(source_functions[[n]]))) &&
    identical(canonical(body(f)),canonical(body(source_functions[[n]])))
},logical(1))]
if (length(bad)) stop("Source/installed function mismatch: ",paste(bad,collapse=", "))
# The on-disk NAMESPACE, not an inherited pkgload namespace, defines exports.
lines <- readLines(file.path(repo,"NAMESPACE"),warn=FALSE)
entries <- lines[startsWith(lines,"export(")]
expected <- substring(entries,8,nchar(entries)-1)
if (!setequal(expected,getNamespaceExports(ns))) stop("Source/installed exports differ")
bundled <- list.files(file.path(repo,"inst"),recursive=TRUE,full.names=TRUE)
bundled <- bundled[grepl("[.](R|py)$",bundled)]
bundled <- bundled[!startsWith(bundled,paste0(file.path(repo,"inst","doc"),"/"))]
for(f in bundled) {
  rel <- substring(f,nchar(file.path(repo,"inst"))+2)
  target <- system.file(rel,package=package,mustWork=TRUE)
  if (!identical(unname(tools::md5sum(f)),unname(tools::md5sum(target)))) stop("Bundled file mismatch: ",rel)
}
sha <- system2("git",c("-C",shQuote(repo),"rev-parse","HEAD"),stdout=TRUE)
dirty <- length(system2("git",c("-C",shQuote(repo),"status","--porcelain"),stdout=TRUE))>0
cat("PASS",package,installed$Version,"source",sha,
    "worktree",if(dirty)"modified" else "clean",
    "functions",length(source_functions),"exports",length(expected),
    "bundled scripts",length(bundled),"installed",find.package(package),"\n")
