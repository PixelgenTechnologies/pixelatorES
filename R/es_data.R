#' Valid diagnostic types for `es_data`
#'
#' @noRd
.es_data_diagnostic_types <- c("pxl_load", "qc_load", "extractor")

#' Create an extractor result with diagnostics
#'
#' Allows an extractor to return a partial value together with non-fatal
#' diagnostics collected while producing it.
#'
#' @param value The extracted slot value.
#' @param diagnostics A list of diagnostics.
#'
#' @return An internal extractor result.
#'
#' @noRd
.new_es_data_extractor_result <- function(value, diagnostics = list()) {
  result <- structure(
    list(value = value, diagnostics = diagnostics),
    class = "es_data_extractor_result"
  )

  return(result)
}

#' Create an `es_data` diagnostic
#'
#' @param type Diagnostic type.
#' @param target Sample alias, pool id, or slot path.
#' @param message Human-readable reason.
#'
#' @return A diagnostic list.
#'
#' @noRd
.new_es_data_diagnostic <- function(type, target, message) {
  diagnostic <- list(
    type = type,
    target = target,
    message = message
  )

  return(diagnostic)
}

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
#'   `sample_aliases`, `effective_samplesheet`, and `file_paths` start as `NULL`
#'   and are set during loading once the samplesheet is known.
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
      sample_aliases = NULL,
      effective_samplesheet = NULL,
      file_paths = NULL,
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
#' nested names under `qc` map to `es_data$qc$...`.
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
    preamble = c("preprocessing.qmd"),
    sections = list(
      list(id = "samples", title = "Samples", child = "samples.qmd"),
      list(
        id = "quality_metrics",
        title = "Quality metrics",
        child = "quality_metrics.qmd"
      ),
      list(
        id = "cell_annotation",
        title = "Cell annotation",
        child = "cell_annotation.qmd"
      ),
      list(id = "abundance", title = "Abundance", child = "abundance.qmd"),
      list(
        id = "spatial",
        title = "Spatial metrics",
        child = "spatial.qmd"
      ),
      list(
        id = "run_settings",
        title = "Run settings",
        child = "run_settings.qmd"
      )
    )
  ))
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

#' Extract merged pixel data
#'
#' Discovers and loads the workflow's pixel files, merges them, and applies
#' report debug-mode downsampling when requested.
#'
#' @param object An `es_data` object.
#'
#' @return A merged Seurat object.
#'
#' @noRd
.extract_pxl_data <- function(object) {
  file_paths <- object$file_paths
  if (is.null(file_paths)) {
    cli_abort("{.field file_paths} must be filled before extracting {.field pxl_data}.")
  }

  pxl_data <- list()
  diagnostics <- list()
  expected_aliases <- object$samplesheet$sample_alias
  available_aliases <- file_paths$data_files$sample_alias[
    file_paths$data_files$sample_alias %in% expected_aliases
  ]
  duplicate_samples <- attr(file_paths, "duplicate_data_samples")
  if (is.null(duplicate_samples)) {
    duplicate_samples <- character()
  }
  aliases_to_load <- c(
    available_aliases,
    setdiff(expected_aliases, c(available_aliases, duplicate_samples))
  )

  for (current_sample_alias in aliases_to_load) {
    sample_file <- file_paths$data_files$filename[
      file_paths$data_files$sample_alias == current_sample_alias
    ]

    if (length(sample_file) == 0) {
      message <- "No PXL file was found."
      cli::cli_warn(
        "Skipping PXL data for sample {.val {current_sample_alias}}: {message}"
      )
      diagnostics <- append(
        diagnostics,
        list(.new_es_data_diagnostic("pxl_load", current_sample_alias, message))
      )
      next
    }

    loaded <- tryCatch(
      list(
        value = ReadPNA_Seurat(
          sample_file,
          load_proximity_scores = FALSE
        ),
        error = NULL
      ),
      error = function(error) {
        return(list(value = NULL, error = error))
      }
    )

    if (!is.null(loaded$error)) {
      message <- conditionMessage(loaded$error)
      cli::cli_warn(
        "Skipping PXL data for sample {.val {current_sample_alias}}: {message}"
      )
      diagnostics <- append(
        diagnostics,
        list(.new_es_data_diagnostic("pxl_load", current_sample_alias, message))
      )
      next
    }

    pxl_data[[current_sample_alias]] <- loaded$value
  }

  if (length(pxl_data) == 0) {
    return(.new_es_data_extractor_result(NULL, diagnostics))
  }

  effective_samplesheet <-
    object$samplesheet %>%
    filter(sample_alias %in% names(pxl_data))

  pxl_data <- merge_data(pxl_data, effective_samplesheet)

  if (isTRUE(object$params$debug_mode)) {
    pxl_data <- downsample_data(
      pxl_data,
      control_markers = object$params$control_markers,
      n_cells = 50,
      n_markers = 21
    )
  }

  return(.new_es_data_extractor_result(pxl_data, diagnostics))
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

#' Extract raw quality control data
#'
#' Discovers and reads per-sample and pool-level QC report files.
#'
#' @param object An `es_data` object.
#'
#' @return A nested list of parsed QC report data.
#'
#' @noRd
.extract_qc_raw <- function(object) {
  file_paths <- object$file_paths
  if (is.null(file_paths)) {
    cli_abort("{.field file_paths} must be filled before extracting {.field qc_raw}.")
  }

  sample_result <- .read_qc_groups_soft(
    files = file_paths$qc_files,
    aliases = sort(object$samplesheet$sample_alias),
    alias_column = "sample_alias",
    target_label = "sample"
  )

  pool_qc_files <- NULL
  pool_diagnostics <- list()
  if ("pool" %in% names(object$samplesheet)) {
    pool_result <- .read_qc_groups_soft(
      files = file_paths$pool_qc_files,
      aliases = sort(unique(object$samplesheet$pool)),
      alias_column = "pool_alias",
      target_label = "pool"
    )
    pool_qc_files <- pool_result$value
    pool_diagnostics <- pool_result$diagnostics
  }

  qc_raw <- list(
    qc_files = sample_result$value,
    pool_qc_files = pool_qc_files
  )
  diagnostics <- c(sample_result$diagnostics, pool_diagnostics)

  return(.new_es_data_extractor_result(qc_raw, diagnostics))
}

#' Read QC file groups without aborting the build
#'
#' Reads all stage files for each expected sample or pool. A missing or
#' unreadable group is skipped and returned as a diagnostic.
#'
#' @param files A data frame of QC file paths.
#' @param aliases Expected sample or pool aliases.
#' @param alias_column Name of the alias column in `files`.
#' @param target_label Label used in warning messages.
#'
#' @return An extractor result containing parsed QC groups and diagnostics.
#'
#' @noRd
.read_qc_groups_soft <- function(
  files,
  aliases,
  alias_column,
  target_label
) {
  values <- list()
  diagnostics <- list()

  for (alias in aliases) {
    if (is.null(files)) {
      alias_files <- NULL
    } else {
      alias_files <- files[files[[alias_column]] == alias, , drop = FALSE]
    }

    if (is.null(alias_files) || nrow(alias_files) == 0) {
      message <- "No QC files were found."
      cli::cli_warn(
        "Skipping QC data for {target_label} {.val {alias}}: {message}"
      )
      diagnostics <- append(
        diagnostics,
        list(.new_es_data_diagnostic("qc_load", alias, message))
      )
      next
    }

    parsed <- tryCatch(
      list(
        value = alias_files$filename %>%
          set_names(alias_files$stage) %>%
          lapply(fload),
        error = NULL
      ),
      error = function(error) {
        return(list(value = NULL, error = error))
      }
    )

    if (!is.null(parsed$error)) {
      message <- conditionMessage(parsed$error)
      cli::cli_warn(
        "Skipping QC data for {target_label} {.val {alias}}: {message}"
      )
      diagnostics <- append(
        diagnostics,
        list(.new_es_data_diagnostic("qc_load", alias, message))
      )
      next
    }

    values[[alias]] <- parsed$value
  }

  return(.new_es_data_extractor_result(values, diagnostics))
}

#' Extract formatted read statistics
#'
#' @param object An `es_data` object.
#'
#' @return A formatted read-statistics table.
#'
#' @noRd
.extract_read_stats <- function(object) {
  metric <- get_read_stats(object$pxl_data)

  return(.format_qc_metric(metric, object$samplesheet))
}

#' Extract formatted sample hashing statistics
#'
#' @param object An `es_data` object.
#'
#' @return Formatted sample hashing statistics, or `NULL`.
#'
#' @noRd
.extract_sample_hash_stats <- function(object) {
  metric <- get_hash_stats(object$pxl_data, object$qc_raw)

  return(.format_qc_metric(metric, object$samplesheet))
}

#' Extract formatted sequencing saturation statistics
#'
#' @param object An `es_data` object.
#'
#' @return A formatted sequencing-saturation table.
#'
#' @noRd
.extract_seq_saturation <- function(object) {
  metric <- get_seq_saturation(object$pxl_data, object$qc_raw)

  return(.format_qc_metric(metric, object$samplesheet))
}

#' Extract formatted crossing-edge statistics
#'
#' @param object An `es_data` object.
#'
#' @return A formatted crossing-edge table.
#'
#' @noRd
.extract_crossing_edges <- function(object) {
  metric <- get_crossing_edges(object$qc_raw)

  return(.format_qc_metric(metric, object$samplesheet))
}

#' Extract formatted degree-distribution statistics
#'
#' @param object An `es_data` object.
#'
#' @return A formatted degree-distribution table.
#'
#' @noRd
.extract_degree_distribution <- function(object) {
  metric <- get_degree_distribution(object$qc_raw)

  return(.format_qc_metric(metric, object$samplesheet))
}

#' Extract formatted denoising statistics
#'
#' @param object An `es_data` object.
#'
#' @return A formatted denoising table, or `NULL`.
#'
#' @noRd
.extract_denoising <- function(object) {
  metric <- get_denoising_data(object$qc_raw)

  return(.format_qc_metric(metric, object$samplesheet))
}

#' Extract formatted detailed denoising statistics
#'
#' @param object An `es_data` object.
#'
#' @return Formatted detailed denoising statistics, or `NULL`.
#'
#' @noRd
.extract_denoising_detail <- function(object) {
  metric <- get_denoising_detail_data(object$pxl_data)

  return(.format_qc_metric(metric, object$samplesheet))
}

#' Extract formatted coreness statistics
#'
#' @param object An `es_data` object.
#'
#' @return A formatted list of coreness tables.
#'
#' @noRd
.extract_coreness <- function(object) {
  metric <- get_coreness_data(object$pxl_data)

  return(.format_qc_metric(metric, object$samplesheet))
}

#' Extract formatted top-marker statistics
#'
#' @param object An `es_data` object.
#'
#' @return A formatted top-marker table.
#'
#' @noRd
.extract_top_markers <- function(object) {
  metric <- get_top_markers(
    object$pxl_data,
    "sample_alias",
    n_markers = c(3, 5)
  )

  return(.format_qc_metric(metric, object$samplesheet))
}

#' Extract processed pixel data
#'
#' @param object An `es_data` object.
#'
#' @return A processed Seurat object.
#'
#' @noRd
.extract_pxl_data_processed <- function(object) {
  return(process_data(object$pxl_data, object$params))
}

#' Derive sample aliases from a samplesheet
#'
#' Creates the named sample and pool alias vector used to order samples in
#' report components.
#'
#' @param samplesheet An Experiment Summary samplesheet.
#'
#' @return A named character vector of sample and pool aliases.
#'
#' @noRd
.sample_aliases_from_samplesheet <- function(samplesheet) {
  sample_aliases <-
    samplesheet %>%
    select(sample, sample_alias) %>%
    deframe()

  if ("pool" %in% names(samplesheet)) {
    pool_aliases <- unique(samplesheet$pool)
    sample_aliases <-
      c(
        sample_aliases,
        set_names(setdiff(pool_aliases, sample_aliases))
      )
  }

  return(sample_aliases)
}

#' Extract filtered proximity scores
#'
#' @param object An `es_data` object.
#'
#' @return A table of filtered proximity scores.
#'
#' @noRd
.extract_proximity <- function(object) {
  proximity <- filter_proximity_scores(
    object$pxl_data_processed,
    object$params,
    sample_levels = object$sample_aliases
  )

  return(proximity)
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

  diagnostic <- .new_es_data_diagnostic(type, target, message)
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
#' An extractor may return an `es_data_extractor_result` to provide a partial
#' value with diagnostics. An unhandled extractor failure records a diagnostic
#' with `type = "extractor"` and leaves its destination slot `NULL` without
#' stopping the remaining extractors. Once `pxl_data_processed` is populated,
#' the raw `pxl_data` slot is cleared to release memory.
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
      error = function(error) {
        return(list(value = NULL, error = error))
      }
    )

    if (!is.null(result$error)) {
      object <- add_es_data_diagnostic(
        object,
        type = "extractor",
        target = target,
        message = conditionMessage(result$error)
      )
    }

    if (inherits(result$value, "es_data_extractor_result")) {
      extractor_result <- result$value
      for (diagnostic in extractor_result$diagnostics) {
        object <- add_es_data_diagnostic(
          object,
          type = diagnostic$type,
          target = diagnostic$target,
          message = diagnostic$message
        )
      }
      result$value <- extractor_result$value
    }

    object <- .set_es_data_slot(object, node_path, result$value)

    if (
      identical(target, "pxl_data_processed") &&
        !is.null(result$value)
    ) {
      object["pxl_data"] <- list(NULL)
    }
  }

  return(object)
}

#' Build Experiment Summary data
#'
#' Creates an `es_data` object, fills the required samplesheet and discovers
#' input file paths once, then runs the remaining workflow extractors.
#' Samplesheet errors are not caught and stop the build; remaining extractor
#' errors are recorded as diagnostics.
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
  object$sample_aliases <-
    .sample_aliases_from_samplesheet(object$samplesheet)
  object <- .fill_es_data_file_paths(object)
  remaining_extractors <- object$extractors[
    setdiff(names(object$extractors), "samplesheet")
  ]

  return(run_es_data_extractors(object, remaining_extractors))
}

#' Discover and store input file paths for an `es_data` object
#'
#' Calls [get_file_paths()] once after the samplesheet is available. Duplicate
#' PXL matches are omitted for the affected sample and recorded as diagnostics.
#'
#' @param object An `es_data` object with `samplesheet` filled.
#'
#' @return The updated `es_data` object.
#'
#' @noRd
.fill_es_data_file_paths <- function(object) {
  file_paths <- get_file_paths(
    data_folder = object$params$data_folder,
    sample_sheet = object$samplesheet,
    on_duplicate_samples = "omit"
  )

  duplicate_samples <- attr(file_paths, "duplicate_data_samples")
  if (is.null(duplicate_samples)) {
    duplicate_samples <- character()
  }
  for (sample_alias in duplicate_samples) {
    message <- "Multiple PXL files matched this sample."
    cli::cli_warn(
      "Skipping PXL data for sample {.val {sample_alias}}: {message}"
    )
    object <- add_es_data_diagnostic(
      object,
      type = "pxl_load",
      target = sample_alias,
      message = message
    )
  }

  object$file_paths <- file_paths

  return(object)
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

#' Summarize an `es_data` slot
#'
#' Creates a compact, human-readable description of the contents of one
#' `es_data` slot for [print.es_data()].
#'
#' @param slot Name of the `es_data` slot.
#' @param value Value stored in the slot.
#'
#' @return A single character string describing `value`.
#'
#' @noRd
.summarize_es_data_slot <- function(slot, value) {
  if (slot == "diagnostics") {
    n_diagnostics <- length(value)
    return(paste(
      n_diagnostics,
      if (n_diagnostics == 1) "diagnostic" else "diagnostics"
    ))
  }

  if (is.null(value) || (is.list(value) && length(value) == 0)) {
    return("not filled")
  }

  if (slot == "params") {
    n_params <- length(value)
    return(paste(n_params, if (n_params == 1) "parameter" else "parameters"))
  }

  if (slot %in% c("samplesheet", "effective_samplesheet")) {
    n_samples <- nrow(value)
    n_columns <- ncol(value)
    return(paste(
      "tibble with",
      n_samples,
      if (n_samples == 1) "sample" else "samples",
      "and",
      n_columns,
      if (n_columns == 1) "column" else "columns"
    ))
  }

  if (slot == "sample_aliases") {
    n_aliases <- length(value)
    return(paste(
      "named character vector with",
      n_aliases,
      if (n_aliases == 1) "alias" else "aliases"
    ))
  }

  if (slot == "file_paths") {
    n_files <- sum(vapply(
      value,
      function(paths) {
        if (is.data.frame(paths)) {
          return(nrow(paths))
        }
        return(0L)
      },
      integer(1)
    ))
    return(paste(
      "list with",
      n_files,
      if (n_files == 1) "discovered file" else "discovered files"
    ))
  }

  if (inherits(value, "Seurat")) {
    return(paste(
      "Seurat object with",
      ncol(value),
      if (ncol(value) == 1) "cell" else "cells",
      "and",
      nrow(value),
      if (nrow(value) == 1) "feature" else "features"
    ))
  }

  if (slot == "qc_raw") {
    n_samples <- length(value$qc_files)
    n_pools <- length(value$pool_qc_files)
    return(paste(
      "raw QC data for",
      n_samples,
      if (n_samples == 1) "sample" else "samples",
      "and",
      n_pools,
      if (n_pools == 1) "pool" else "pools"
    ))
  }

  if (slot == "qc") {
    n_metrics <- sum(!vapply(value, is.null, logical(1)))
    return(paste(
      "list with",
      n_metrics,
      if (n_metrics == 1) "formatted metric" else "formatted metrics"
    ))
  }

  if (slot == "proximity" && is.data.frame(value)) {
    n_scores <- nrow(value)
    return(paste(
      "tibble with",
      n_scores,
      if (n_scores == 1) "proximity score" else "proximity scores"
    ))
  }

  if (slot == "extractors") {
    n_extractors <- sum(vapply(
      unlist(value, recursive = TRUE),
      is.function,
      logical(1)
    ))
    return(paste(
      "workflow registry with",
      n_extractors,
      if (n_extractors == 1) "extractor" else "extractors"
    ))
  }

  return(paste(class(value), collapse = "/"))
}

#' Print an `es_data` object
#'
#' Prints the workflow and a concise description of the contents of every
#' `es_data` slot.
#'
#' @param x An `es_data` object.
#' @param ... Not used.
#'
#' @return `x`, invisibly.
#'
#' @export
print.es_data <- function(x, ...) {
  cli::cli_h1("es_data")
  cli::cli_text("Workflow: {.val {x$params$workflow}}")
  cli::cli_ul()
  for (slot in names(x)) {
    summary <- .summarize_es_data_slot(slot, x[[slot]])
    if (summary == "not filled") {
      cli::cli_li("{.field {slot}}: {.emph {summary}}")
    } else {
      cli::cli_li("{.field {slot}}: {summary}")
    }
  }
  cli::cli_end()

  return(invisible(x))
}
