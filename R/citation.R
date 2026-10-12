# Part of the InfoMergeCluster package.
# Author: Komlan Midodzi Noukpoape <noukpoape@gmail.com>

# BibTeX entry of the package (Google Scholar format); also installed as
# system.file("CITATION.bib", package = "InfoMergeCluster")
.im_bibtex <- paste(
  "@manual{noukpoape2026infomergecluster,",
  "  title={InfoMergeCluster: Clustering of Sparse Longitudinal Data with Automatic Selection of the Number of Clusters},",
  "  author={Noukpoape, Komlan Midodzi},",
  "  year={2026},",
  "  note={R package version 1.0.0},",
  "  url={https://github.com/noukpoape/InfoMergeCluster}",
  "}",
  sep = "\n"
)

.im_cite_text <- paste(
  "Noukpoape, K. M. (2026). InfoMergeCluster: Clustering of Sparse Longitudinal",
  "Data with Automatic Selection of the Number of Clusters. R package version 1.0.0.",
  "https://github.com/noukpoape/InfoMergeCluster"
)

#' How to cite InfoMergeCluster
#'
#' Prints the reference to cite when using the package, either as a BibTeX
#' entry in Google Scholar format or as plain text.
#'
#' The BibTeX entry is also installed with the package, in the file
#' \code{system.file("CITATION.bib", package = "InfoMergeCluster")}, and
#' \code{citation("InfoMergeCluster")} returns the same reference.
#'
#' @param format \code{"bibtex"} (default) or \code{"text"}.
#' @param file optional path of a \code{.bib} file to which the BibTeX entry
#'   is appended (only with \code{format = "bibtex"}).
#' @return The reference as a character string, invisibly.
#' @examples
#' infomergecluster_citation()
#' infomergecluster_citation("text")
#' @export
infomergecluster_citation <- function(format = c("bibtex", "text"), file = NULL) {
  format <- match.arg(format)
  ref <- if (format == "bibtex") .im_bibtex else .im_cite_text
  cat(ref, "\n", sep = "")
  if (!is.null(file)) {
    if (format != "bibtex") stop("'file' can only be used with format = \"bibtex\".")
    cat(ref, "\n\n", file = file, sep = "", append = TRUE)
    message(sprintf("BibTeX entry appended to '%s'.", file))
  }
  invisible(ref)
}

.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
    "InfoMergeCluster 1.0.0. To cite the package in publications, ",
    "use infomergecluster_citation()."
  )
}
