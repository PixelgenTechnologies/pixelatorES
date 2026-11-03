#' Extract formatted sequencing saturation statistics for edgelist_aggregate
#'
#' @param object An `es_data` object.
#'
#' @return A formatted sequencing-saturation table.
#'
#' @noRd
.extract_edgelist_aggregate_seq_saturation <- function(object) {
  metric <- get_seq_saturation(
    object$pxl_data,
    object$qc_raw,
    read_stage = "edgelist",
    collapse_stage = "aggregate"
  )

  return(.format_qc_metric(metric, object$samplesheet))
}

#' Extract formatted degree-distribution statistics for edgelist_aggregate
#'
#' @param object An `es_data` object.
#'
#' @return A formatted degree-distribution table.
#'
#' @noRd
.extract_edgelist_aggregate_degree_distribution <- function(object) {
  metric <- get_degree_distribution(object$qc_raw, stage = "aggregate")

  return(.format_qc_metric(metric, object$samplesheet))
}

#' Extractor registry for the edgelist_aggregate workflow
#'
#' Nested named list of functions. Top-level names map to `es_data` slots;
#' nested names under `qc` map to `es_data$qc$...`. The generic extractors
#' live in [es_data.R](R/es_data.R); the two that read `edgelist` and
#' `aggregate` call the parameterised `get_*()` functions with those stage
#' names.
#'
#' @return A nested named list of functions.
#'
#' @noRd
.edgelist_aggregate_extractors <- function() {
  extractors <- list(
    samplesheet = .extract_samplesheet,
    pxl_data = .extract_pxl_data,
    effective_samplesheet = .extract_effective_samplesheet,
    qc_raw = .extract_qc_raw,
    qc = list(
      read_stats = .extract_read_stats,
      sample_hash_stats = .extract_sample_hash_stats,
      seq_saturation = .extract_edgelist_aggregate_seq_saturation,
      crossing_edges = .extract_crossing_edges,
      degree_distribution = .extract_edgelist_aggregate_degree_distribution,
      denoising = .extract_denoising,
      denoising_detail = .extract_denoising_detail,
      coreness = .extract_coreness,
      top_markers = .extract_top_markers,
      graph_collapse = .extract_graph_collapse
    ),
    pxl_data_processed = .extract_pxl_data_processed,
    proximity = .extract_proximity
  )

  return(extractors)
}

#' Report recipe for the edgelist_aggregate workflow
#'
#' Paths are relative to `inst/quarto/`.
#'
#' @return A report recipe list.
#'
#' @noRd
.edgelist_aggregate_report <- function() {
  return(list(
    preamble = c("shared/preprocessing.qmd"),
    sections = list(
      list(id = "samples", title = "Samples", child = "shared/samples.qmd"),
      list(
        id = "quality_metrics",
        title = "Quality metrics",
        child = "workflows/edgelist_aggregate/quality_metrics.qmd"
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

#' Pipeline stage vocabulary for the edgelist_aggregate workflow
#'
#' Stage names emitted by an nf-core/pixelator run whose read-processing front
#' end is `edgelist` + `aggregate` rather than `amplicon` + `demux` +
#' `collapse`. `edgelist` reports the merged amplicon and demultiplexing
#' statistics, and `aggregate` reports the read-to-molecule collapse.
#' `graph-collapse` runs on the `graph` output and merges nearby nodes with
#' near-identical sequences; its reports are staged under `graph_collapse`.
#'
#' `denoise` is pool-level here because it runs before `sample_calling`, whereas
#' in `amplicon_demux` it runs after and is therefore per-sample.
#'
#' @return A stage vocabulary list.
#'
#' @noRd
.edgelist_aggregate_stages <- function() {
  return(list(
    all = c(
      "edgelist", "aggregate",
      "graph", "graph_collapse", "denoise",
      "sample_calling",
      "analysis", "post_analysis", "layout"
    ),
    pool = c(
      "edgelist", "aggregate",
      "graph", "graph_collapse", "denoise"
    ),
    pxl_preference = c("graph", "analysis", "post_analysis", "layout")
  ))
}
