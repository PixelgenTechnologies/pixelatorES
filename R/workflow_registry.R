#' Registered Experiment Summary workflows
#'
#' Mutable package-local registry mapping workflow identifiers to workflow
#' definitions. Each definition contains zero-argument factory functions for
#' extractors, the report recipe, and the pipeline stage vocabulary.
#'
#' @noRd
.es_data_workflow_registry <- new.env(parent = emptyenv())

#' Required elements of a workflow stage vocabulary
#'
#' @noRd
.es_workflow_stage_fields <- c("all", "pool", "pxl_preference")

#' Resolve a report child path for existence checks
#'
#' Absolute paths are returned unchanged. Relative paths are resolved against
#' `inst/quarto/` in pixelatorES (built-in workflow convention).
#'
#' @param path A child document path.
#'
#' @return The path to check with [file.exists()].
#'
#' @noRd
.resolve_es_workflow_report_path <- function(path) {
  if (grepl("^(~|/|[A-Za-z]:[/\\\\]|\\\\\\\\)", path)) {
    return(path)
  }

  quarto_root <- system.file("quarto", package = "pixelatorES")
  if (!nzchar(quarto_root)) {
    cli_abort("Could not locate installed {.pkg pixelatorES} Quarto directory.")
  }

  return(file.path(quarto_root, path))
}

#' Validate that report recipe child paths exist
#'
#' @param report A structurally valid report recipe.
#'
#' @return `report`, invisibly.
#'
#' @noRd
.validate_es_workflow_report_paths <- function(report) {
  paths <- c(
    report$preamble,
    vapply(report$sections, function(section) {
      return(section$child)
    }, character(1))
  )
  resolved <- vapply(paths, .resolve_es_workflow_report_path, character(1))
  missing <- paths[!file.exists(resolved)]
  if (length(missing) > 0) {
    cli_abort(c(
      "Report recipe references missing Quarto child documents.",
      "x" = "Missing: {.val {unique(missing)}}.",
      "i" = paste(
        "Built-in workflows use paths relative to {.path inst/quarto/};",
        "extension packages should register absolute paths from {.fn system.file}."
      )
    ))
  }

  return(invisible(report))
}

#' Validate an Experiment Summary report recipe
#'
#' @param report A report recipe list.
#'
#' @return The validated `report`.
#'
#' @noRd
.validate_es_workflow_report <- function(report) {
  pixelatorR:::assert_class(report, "list")
  if (is.null(names(report)) || any(!nzchar(names(report)))) {
    cli_abort("{.arg report} must be a named list.")
  }

  unexpected <- setdiff(names(report), c("preamble", "sections"))
  if (length(unexpected) > 0) {
    cli_abort(c(
      "{.arg report} contains unexpected elements.",
      "x" = "Unexpected: {.val {unexpected}}.",
      "i" = "Allowed elements: {.val preamble}, {.val sections}."
    ))
  }

  if (is.null(report$preamble)) {
    cli_abort("{.arg report} must contain {.field preamble}.")
  }
  # assert_vector(..., n) requires at least n elements (not exactly n).
  pixelatorR:::assert_vector(report$preamble, "character", n = 1)
  if (any(!nzchar(report$preamble))) {
    cli_abort("{.arg report$preamble} must not contain empty paths.")
  }

  if (is.null(report$sections)) {
    cli_abort("{.arg report} must contain {.field sections}.")
  }
  pixelatorR:::assert_class(report$sections, "list")
  if (length(report$sections) == 0) {
    cli_abort("{.arg report$sections} must contain at least one section.")
  }

  section_ids <- character(length(report$sections))
  for (i in seq_along(report$sections)) {
    section <- report$sections[[i]]
    pixelatorR:::assert_class(section, "list")
    for (field in c("id", "title", "child")) {
      if (is.null(section[[field]])) {
        cli_abort(c(
          "Report section {.val {i}} is missing {.field {field}}.",
          "i" = "Each section needs {.field id}, {.field title}, and {.field child}."
        ))
      }
      pixelatorR:::assert_single_value(section[[field]], "string")
      if (!nzchar(section[[field]])) {
        cli_abort("Report section {.val {i}} has an empty {.field {field}}.")
      }
    }
    section_ids[[i]] <- section$id
  }

  if (anyDuplicated(section_ids) > 0) {
    duplicated_ids <- unique(section_ids[duplicated(section_ids)])
    cli_abort(c(
      "{.arg report$sections} contains duplicated section ids.",
      "x" = "Duplicated: {.val {duplicated_ids}}."
    ))
  }

  .validate_es_workflow_report_paths(report)

  return(report)
}

#' Validate a workflow report factory
#'
#' @param report A zero-argument function returning a report recipe.
#'
#' @return `report`, invisibly.
#'
#' @noRd
.validate_es_workflow_report_factory <- function(report) {
  pixelatorR:::assert_class(report, "function")
  .validate_es_workflow_report(report())

  return(invisible(report))
}

#' Validate a workflow extractor factory
#'
#' @param extractors A zero-argument function returning a nested named list of
#'   extractor functions.
#'
#' @return `extractors`, invisibly.
#'
#' @noRd
.validate_es_workflow_extractors <- function(extractors) {
  pixelatorR:::assert_class(extractors, "function")

  extractor_list <- extractors()
  pixelatorR:::assert_class(extractor_list, "list")

  return(invisible(extractors))
}

#' Validate an Experiment Summary workflow stage vocabulary
#'
#' @param stages A stage vocabulary list.
#'
#' @return The validated `stages`.
#'
#' @noRd
.validate_es_workflow_stages <- function(stages) {
  pixelatorR:::assert_class(stages, "list")
  if (is.null(names(stages)) || any(!nzchar(names(stages)))) {
    cli_abort("{.arg stages} must be a named list.")
  }

  missing <- setdiff(.es_workflow_stage_fields, names(stages))
  if (length(missing) > 0) {
    cli_abort(c(
      "{.arg stages} is missing required elements.",
      "x" = "Missing: {.val {missing}}.",
      "i" = "Required elements: {.val {.es_workflow_stage_fields}}."
    ))
  }

  unexpected <- setdiff(names(stages), .es_workflow_stage_fields)
  if (length(unexpected) > 0) {
    cli_abort(c(
      "{.arg stages} contains unexpected elements.",
      "x" = "Unexpected: {.val {unexpected}}.",
      "i" = "Allowed elements: {.val {.es_workflow_stage_fields}}."
    ))
  }

  # assert_vector(..., n) requires at least n elements (not exactly n), so
  # n = 0 admits an empty `pool` for pipelines without pool-level QC.
  pixelatorR:::assert_vector(stages$all, "character", n = 1)
  pixelatorR:::assert_vector(stages$pool, "character", n = 0)
  pixelatorR:::assert_vector(stages$pxl_preference, "character", n = 1)

  for (field in .es_workflow_stage_fields) {
    field_label <- paste0("stages$", field)
    if (any(!nzchar(stages[[field]]))) {
      cli_abort("{.arg {field_label}} must not contain empty stage names.")
    }
    if (anyDuplicated(stages[[field]]) > 0) {
      duplicated_stages <- unique(stages[[field]][duplicated(stages[[field]])])
      cli_abort(c(
        "{.arg {field_label}} contains duplicated stages.",
        "x" = "Duplicated: {.val {duplicated_stages}}."
      ))
    }
  }

  for (field in c("pool", "pxl_preference")) {
    field_label <- paste0("stages$", field)
    unknown <- setdiff(stages[[field]], stages$all)
    if (length(unknown) > 0) {
      cli_abort(c(
        "{.arg {field_label}} contains stages missing from {.field all}.",
        "x" = "Unknown: {.val {unknown}}."
      ))
    }
  }

  return(stages)
}

#' Validate a workflow stage factory
#'
#' @param stages A zero-argument function returning a stage vocabulary.
#'
#' @return `stages`, invisibly.
#'
#' @noRd
.validate_es_workflow_stages_factory <- function(stages) {
  pixelatorR:::assert_class(stages, "function")
  .validate_es_workflow_stages(stages())

  return(invisible(stages))
}

#' Register an Experiment Summary workflow
#'
#' Registers a workflow identifier with its extractor factory, Quarto report
#' recipe, and pipeline stage vocabulary. Extension packages should call this
#' from their `.onLoad()` hook. Workflows are the only supported way to build
#' Experiment Summary data and render the report.
#'
#' @param name A unique workflow identifier.
#' @param extractors A zero-argument function returning a nested named list of
#'   extractor functions. The factory is called at registration time to verify
#'   that it returns a list.
#' @param report A zero-argument function returning the report recipe used by
#'   the Quarto shell. The recipe is a named list with:
#'   - `preamble`: non-empty character vector of child document paths knitted
#'     before the tabset (for example data ingestion).
#'   - `sections`: non-empty list of tabset sections, each a list with `id`,
#'     `title`, and `child`.
#'   Built-in workflows use paths relative to `inst/quarto/`. Extension packages
#'   should register absolute paths from `system.file()`. All referenced child
#'   paths must exist at registration time.
#' @param stages A zero-argument function returning the pipeline stage
#'   vocabulary used for input file discovery. The vocabulary is a named list
#'   with:
#'   - `all`: non-empty character vector of every stage the pipeline can emit.
#'   - `pool`: character vector of stages whose QC files are pool-level. May be
#'     empty for pipelines without pool-level QC.
#'   - `pxl_preference`: non-empty character vector of stages ordered by
#'     preference when several stages produce a PXL file for a sample.
#'   `pool` and `pxl_preference` must be subsets of `all`.
#' @param overwrite If `TRUE`, replace an existing registration for `name`.
#'   Defaults to `FALSE`.
#'
#' @return `name`, invisibly.
#'
#' @export
register_es_data_workflow <- function(
  name,
  extractors,
  report,
  stages,
  overwrite = FALSE
) {
  pixelatorR:::assert_single_value(name, "string")
  pixelatorR:::assert_single_value(overwrite, "bool")
  .validate_es_workflow_extractors(extractors)
  .validate_es_workflow_report_factory(report)
  .validate_es_workflow_stages_factory(stages)

  if (
    !overwrite &&
      exists(name, envir = .es_data_workflow_registry, inherits = FALSE)
  ) {
    cli_abort("Workflow {.val {name}} is already registered.")
  }

  assign(
    name,
    list(extractors = extractors, report = report, stages = stages),
    envir = .es_data_workflow_registry
  )

  return(invisible(name))
}

#' List registered Experiment Summary workflows
#'
#' @return A character vector of registered workflow identifiers.
#'
#' @export
list_es_data_workflows <- function() {
  return(ls(envir = .es_data_workflow_registry, all.names = TRUE, sorted = TRUE))
}

#' Get a registered workflow definition
#'
#' @param name A workflow identifier.
#'
#' @return A list with `extractors`, `report`, and `stages`.
#'
#' @noRd
.get_es_workflow_definition <- function(name) {
  pixelatorR:::assert_single_value(name, "string")

  if (!exists(name, envir = .es_data_workflow_registry, inherits = FALSE)) {
    cli_abort(c(
      "Unsupported {.arg workflow}: {.val {name}}.",
      "i" = "Registered workflows: {.val {list_es_data_workflows()}}."
    ))
  }

  return(get(
    name,
    envir = .es_data_workflow_registry,
    inherits = FALSE
  ))
}

#' Get extractors for a registered workflow
#'
#' @param name A workflow identifier.
#'
#' @return A nested named list of extractor functions.
#'
#' @noRd
.get_es_data_extractors <- function(name) {
  return(.get_es_workflow_definition(name)$extractors())
}

#' Get the report recipe for a registered workflow
#'
#' Returns the Quarto report recipe registered for `name`. Built-in workflows
#' use paths relative to `inst/quarto/`; extension packages may register
#' absolute paths.
#'
#' @param name A workflow identifier.
#'
#' @return A report recipe list with `preamble` and `sections`.
#'
#' @export
get_es_workflow_report <- function(name) {
  return(.get_es_workflow_definition(name)$report())
}

#' Get the pipeline stage vocabulary for a registered workflow
#'
#' Returns the stage vocabulary registered for `name`. It drives input file
#' discovery: which stage names are legal, which stages carry pool-level QC,
#' and how PXL files are preferred when several stages produce one.
#'
#' @param name A workflow identifier.
#'
#' @return A stage vocabulary list with `all`, `pool`, and `pxl_preference`.
#'
#' @export
get_es_workflow_stages <- function(name) {
  return(.get_es_workflow_definition(name)$stages())
}
