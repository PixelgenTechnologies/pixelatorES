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
#' @param params A list of Experiment Summary parameters. If `workflow` is
#'   absent, it defaults to `"amplicon_demux"`.
#'
#' @return An `es_data` object with extractors attached but data slots unfilled.
#'   `effective_samplesheet` starts as `NULL` and is set during loading once
#'   successful samples are known.
#'
#' @keywords internal
new_es_data <- function(params) {
  pixelatorR:::assert_class(params, c("list", "knit_param_list"))

  workflow <- params$workflow
  if (is.null(workflow)) {
    workflow <- "amplicon_demux"
    params$workflow <- workflow
  }
  pixelatorR:::assert_single_value(workflow, "string")

  object <- structure(
    list(
      params = params,
      diagnostics = list(),
      samplesheet = NULL,
      effective_samplesheet = NULL,
      pxl_data = NULL,
      qc_raw = NULL,
      qc = list(),
      pxl_data_processed = NULL,
      proximity = NULL,
      extractors = .get_es_data_extractors(workflow)
    ),
    class = c("es_data", "list")
  )

  return(object)
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
  extractors <- list(
    samplesheet = .extract_samplesheet,
    pxl_data = function(object) NULL,
    effective_samplesheet = .extract_effective_samplesheet,
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

  return(extractors)
}

#' Extract the experiment samplesheet
#'
#' @param object An `es_data` object containing `params$sample_sheet`.
#'
#' @return The parsed samplesheet.
#'
#' @noRd
.extract_samplesheet <- function(object) {
  return(read_samplesheet(object$params$sample_sheet))
}

#' Extract the effective samplesheet
#'
#' Filters the original samplesheet to samples represented in `pxl_data`. If
#' no pxl data loaded, returns an empty samplesheet with the original columns.
#'
#' @param object An `es_data` object.
#'
#' @return A samplesheet containing samples represented in `pxl_data`.
#'
#' @noRd
.extract_effective_samplesheet <- function(object) {
  loaded_samples <-
    if (is.null(object$pxl_data)) {
      character()
    } else {
      FetchData(object$pxl_data, "sample_alias") %>%
        pull(sample_alias) %>%
        as.character() %>%
        unique()
    }

  effective_samplesheet <-
    object$samplesheet %>%
    filter(sample_alias %in% loaded_samples)

  return(effective_samplesheet)
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

  return(object)
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

  return(object)
}

#' Build Experiment Summary data
#'
#' Creates an `es_data` object, fills the required samplesheet, then runs the
#' remaining workflow extractors. Samplesheet errors are not caught and stop
#' the build; remaining extractor errors are recorded as diagnostics.
#'
#' @param params A list of Experiment Summary parameters.
#'
#' @return An `es_data` object.
#'
#' @export
build_es_data <- function(params) {
  object <- new_es_data(params)

  if (is.null(object$extractors$samplesheet)) {
    cli_abort(
      "Workflow {.val {object$params$workflow}} must define a samplesheet extractor."
    )
  }

  object$samplesheet <- object$extractors$samplesheet(object)
  remaining_extractors <- object$extractors[
    setdiff(names(object$extractors), "samplesheet")
  ]

  return(run_es_data_extractors(object, remaining_extractors))
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

  return(object)
}

#' Print an `es_data` object
#'
#' Prints the workflow, how many data slots are populated, and the number of
#' recorded diagnostics.
#'
#' @param x An `es_data` object.
#' @param ... Not used.
#'
#' @return `x`, invisibly.
#'
#' @export
print.es_data <- function(x, ...) {
  data_slots <- setdiff(names(x), c("params", "diagnostics", "extractors"))
  populated <- vapply(x[data_slots], Negate(is.null), logical(1))

  cat("<es_data:", x$params$workflow, ">\n", sep = "")
  cat("Populated slots:", sum(populated), "of", length(populated), "\n")
  cat("Diagnostics:", length(x$diagnostics), "\n")
  invisible(x)
}
