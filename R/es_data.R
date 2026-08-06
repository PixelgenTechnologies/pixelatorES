#' Supported Experiment Summary workflows
#'
#' @noRd
.es_data_workflows <- c("amplicon_demux")

#' Valid diagnostic types for `es_data`
#'
#' @noRd
.es_data_diagnostic_types <- c("pxl_load", "qc_load", "extractor")

#' Create an Experiment Summary data object
#'
#' `es_data` is the single data object passed through Experiment Summary
#' preprocessing. Its slots are filled by workflow-specific extractors attached
#' at construction time.
#'
#' @param workflow The workflow used to select extractors. Defaults to
#'   `"amplicon_demux"`.
#' @param samplesheet The original experiment samplesheet, or `NULL` until set.
#' @param meta A named list of build metadata.
#'
#' @return An `es_data` object with extractors attached but data slots unfilled.
#'   `effective_samplesheet` starts as `NULL` and is set during loading once
#'   successful samples are known.
#'
#' @keywords internal
new_es_data <- function(
  workflow = "amplicon_demux",
  samplesheet = NULL,
  meta = list()
) {
  pixelatorR:::assert_single_value(workflow, "string")
  pixelatorR:::assert_x_in_y(workflow, .es_data_workflows)
  pixelatorR:::assert_class(samplesheet, "data.frame", allow_null = TRUE)
  pixelatorR:::assert_class(meta, "list")

  meta$workflow <- workflow

  structure(
    list(
      meta = meta,
      diagnostics = list(),
      samplesheet = samplesheet,
      effective_samplesheet = NULL,
      pxl_data = NULL,
      qc_raw = NULL,
      qc = list(),
      pxl_data_processed = NULL,
      proximity = NULL,
      extractors = .es_data_extractors_for_workflow(workflow)
    ),
    class = c("es_data", "list")
  )
}

#' Select extractors for a workflow
#'
#' @param workflow Workflow name.
#'
#' @return A nested named list of extractor functions.
#'
#' @noRd
.es_data_extractors_for_workflow <- function(workflow) {
  switch(
    workflow,
    amplicon_demux = .amplicon_demux_extractors(),
    cli_abort(c(
      "Unsupported {.arg workflow}: {.val {workflow}}.",
      "i" = "Supported workflows: {.val {.es_data_workflows}}."
    ))
  )
}

#' Extractor registry for the amplicon_demux workflow
#'
#' Nested named list of functions. Top-level names map to `es_data` slots;
#' nested names under `qc` map to `es_data$qc$...`. Phase 1b replaces stubs
#' with adapters around existing loaders and getters.
#'
#' @return A nested named list of functions.
#'
#' @noRd
.amplicon_demux_extractors <- function() {
  list(
    pxl_data = function(object) NULL,
    qc_raw = function(object) NULL,
    qc = list(
      read_stats = function(object) NULL,
      sample_hash_stats = function(object) NULL,
      seq_saturation = function(object) NULL,
      crossing_edges = function(object) NULL,
      degree_distribution = function(object) NULL,
      denoising = function(object) NULL,
      denoising_detail = function(object) NULL,
      coreness = function(object) NULL,
      top_markers = function(object) NULL
    ),
    pxl_data_processed = function(object) NULL,
    proximity = function(object) NULL
  )
}

#' Add a diagnostic to an `es_data` object
#'
#' Soft-fail issues are recorded as diagnostics. Membership in
#' `object$diagnostics` implies a soft failure.
#'
#' @param object An `es_data` object.
#' @param type Diagnostic type: `"pxl_load"`, `"qc_load"`, or `"extractor"`.
#' @param target Sample alias, pool id, or slot path that failed
#'   (e.g. `"proximity"` or `"qc$crossing_edges"`).
#' @param message Human-readable reason.
#'
#' @return The updated `es_data` object.
#'
#' @keywords internal
add_es_data_diagnostic <- function(object, type, target, message) {
  pixelatorR:::assert_class(object, "es_data")
  pixelatorR:::assert_single_value(type, "string")
  pixelatorR:::assert_is_one_of(type, .es_data_diagnostic_types)
  pixelatorR:::assert_single_value(target, "string")
  pixelatorR:::assert_single_value(message, "string")

  diagnostic <- list(
    type = type,
    target = target,
    message = message
  )
  object$diagnostics <- append(object$diagnostics, list(diagnostic))
  object
}

#' Run registered `es_data` extractors
#'
#' Walks a nested named list of extractor functions depth-first. Each leaf must
#' be a function that takes the current `es_data` and returns the value for that
#' slot. Nested lists (e.g. under `qc`) fill nested slots such as
#' `object$qc$read_stats`.
#'
#' An extractor failure records a diagnostic with `type = "extractor"` and
#' leaves its destination slot `NULL` without stopping the remaining
#' extractors.
#'
#' @param object An `es_data` object.
#' @param extractors Nested named list of functions. Defaults to
#'   `object$extractors`. Used recursively for nested registries.
#' @param path Character vector of slot names from the root. Used recursively
#'   to build diagnostic targets such as `"qc$crossing_edges"`.
#'
#' @return The updated `es_data` object.
#'
#' @keywords internal
run_es_data_extractors <- function(
  object,
  extractors = object$extractors,
  path = character()
) {
  pixelatorR:::assert_class(object, "es_data")

  for (name in names(extractors)) {
    node <- extractors[[name]]
    node_path <- c(path, name)
    target <- paste(node_path, collapse = "$")

    if (is.list(node) && !is.function(node)) {
      object <- run_es_data_extractors(object, node, node_path)
      next
    }

    result <- tryCatch(
      list(value = node(object), error = NULL),
      error = function(error) list(value = NULL, error = error)
    )

    if (!is.null(result$error)) {
      object <- add_es_data_diagnostic(
        object,
        type = "extractor",
        target = target,
        message = conditionMessage(result$error)
      )
    }

    object <- .set_es_data_slot(object, node_path, result$value)
  }

  object
}

#' Build Experiment Summary data
#'
#' Creates an `es_data` object for `workflow` and runs its extractors. The
#' workflow selects which extractor registry fills the bag.
#'
#' @param workflow The workflow to build. Defaults to `"amplicon_demux"`.
#' @param samplesheet The original experiment samplesheet, or `NULL`.
#' @param meta A named list of build metadata.
#'
#' @return An `es_data` object.
#'
#' @export
build_es_data <- function(
  workflow = "amplicon_demux",
  samplesheet = NULL,
  meta = list()
) {
  object <- new_es_data(
    workflow = workflow,
    samplesheet = samplesheet,
    meta = meta
  )
  run_es_data_extractors(object)
}

#' Set a (possibly nested) slot on an `es_data`-like list
#'
#' @param object A list or `es_data` object.
#' @param slot Character vector of nested names.
#' @param value Value to store.
#'
#' @return The updated object.
#'
#' @noRd
.set_es_data_slot <- function(object, slot, value) {
  slot_name <- slot[[1]]

  if (length(slot) == 1) {
    object[slot_name] <- list(value)
    return(object)
  }

  child <- object[[slot_name]]
  if (is.null(child)) {
    child <- list()
  }
  if (!is.list(child)) {
    cli_abort(c(
      "Cannot fill slot {.field {paste(slot, collapse = '$')}}.",
      "x" = "{.field {slot_name}} is not a list."
    ))
  }

  object[slot_name] <- list(.set_es_data_slot(child, slot[-1], value))
  object
}

#' @export
print.es_data <- function(x, ...) {
  data_slots <- setdiff(names(x), c("meta", "diagnostics", "extractors"))
  populated <- vapply(x[data_slots], Negate(is.null), logical(1))

  cat("<es_data:", x$meta$workflow, ">\n", sep = "")
  cat("Populated slots:", sum(populated), "of", length(populated), "\n")
  cat("Diagnostics:", length(x$diagnostics), "\n")
  invisible(x)
}
