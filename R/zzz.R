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
    stages = .amplicon_demux_stages,
    overwrite = TRUE
  )

  reg.finalizer(
    .psock_pool,
    function(e) {
      if (!is.null(e$cl)) {
        try(parallel::stopCluster(e$cl), silent = TRUE)
      }
    },
    onexit = TRUE
  )

  return(invisible(NULL))
}
