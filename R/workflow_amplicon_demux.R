#' Extractor registry for the amplicon_demux workflow
#'
#' Nested named list of functions. Top-level names map to `es_data` slots;
#' nested names under `qc` map to `es_data$qc$...`. The extractor
#' implementations live in [es_data.R](R/es_data.R).
#'
#' @return A nested named list of functions.
#'
#' @noRd
.amplicon_demux_extractors <- function() {
  extractors <- list(
    samplesheet = .extract_samplesheet,
    pxl_data = .extract_pxl_data,
    effective_samplesheet = .extract_effective_samplesheet,
    qc_raw = .extract_qc_raw,
    qc = list(
      read_stats = .extract_read_stats,
      sample_hash_stats = .extract_sample_hash_stats,
      seq_saturation = .extract_seq_saturation,
      crossing_edges = .extract_crossing_edges,
      degree_distribution = .extract_degree_distribution,
      denoising = .extract_denoising,
      denoising_detail = .extract_denoising_detail,
      coreness = .extract_coreness,
      top_markers = .extract_top_markers
    ),
    pxl_data_processed = .extract_pxl_data_processed,
    proximity = .extract_proximity
  )

  return(extractors)
}

#' Report recipe for the amplicon_demux workflow
#'
#' Paths are relative to `inst/quarto/`.
#'
#' @return A report recipe list.
#'
#' @noRd
.amplicon_demux_report <- function() {
  return(list(
    preamble = c("shared/preprocessing.qmd"),
    sections = list(
      list(id = "samples", title = "Samples", child = "shared/samples.qmd"),
      list(
        id = "quality_metrics",
        title = "Quality metrics",
        child = "workflows/amplicon_demux/quality_metrics.qmd"
      ),
      list(
        id = "cell_annotation",
        title = "Cell annotation",
        child = "workflows/amplicon_demux/cell_annotation.qmd"
      ),
      list(
        id = "abundance",
        title = "Abundance",
        child = "workflows/amplicon_demux/abundance.qmd"
      ),
      list(
        id = "spatial",
        title = "Spatial metrics",
        child = "workflows/amplicon_demux/spatial.qmd"
      ),
      list(
        id = "run_info",
        title = "Run info",
        child = "shared/run_info.qmd"
      )
    )
  ))
}

#' Pipeline stage vocabulary for the amplicon_demux workflow
#'
#' Stage names emitted by nf-core/pixelator, the subset whose QC files are
#' pool-level, and the order in which stages are preferred when several
#' produce a PXL file for the same sample.
#'
#' @return A stage vocabulary list.
#'
#' @noRd
.amplicon_demux_stages <- function() {
  return(list(
    all = c(
      "amplicon", "collapse",
      "demux", "denoise",
      "graph", "analysis",
      "sample_calling",
      "post_analysis", "layout"
    ),
    pool = c(
      "amplicon", "collapse",
      "demux", "graph"
    ),
    pxl_preference = c("graph", "analysis", "post_analysis", "layout")
  ))
}
