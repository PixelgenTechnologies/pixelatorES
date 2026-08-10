#' Register built-in workflows when pixelatorES loads
#'
#' @param libname Path to the package library.
#' @param pkgname Package name.
#'
#' @return `NULL`, invisibly.
#'
#' @noRd
.onLoad <- function(libname, pkgname) {
  register_es_data_workflow(
    name = "amplicon_demux",
    extractors = .amplicon_demux_extractors,
    report = .amplicon_demux_report,
    overwrite = TRUE
  )

  return(invisible(NULL))
}
