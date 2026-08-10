#' Registered Experiment Summary workflows
#'
#' Mutable package-local registry mapping workflow identifiers to workflow
#' definitions. Each definition contains zero-argument factory functions for
#' extractors and the report recipe.
#'
#' @noRd
.es_data_workflow_registry <- new.env(parent = emptyenv())

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

#' Register an Experiment Summary workflow
#'
#' Registers a workflow identifier with its extractor factory and Quarto report
#' recipe. Extension packages should call this from their `.onLoad()` hook.
#' Workflows are the only supported way to build Experiment Summary data and
#' render the report.
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
#'   should register absolute paths from `system.file()`.
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
  overwrite = FALSE
) {
  pixelatorR:::assert_single_value(name, "string")
  pixelatorR:::assert_single_value(overwrite, "bool")
  .validate_es_workflow_extractors(extractors)
  .validate_es_workflow_report_factory(report)

  if (
    !overwrite &&
      exists(name, envir = .es_data_workflow_registry, inherits = FALSE)
  ) {
    cli_abort("Workflow {.val {name}} is already registered.")
  }

  assign(
    name,
    list(extractors = extractors, report = report),
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
#' @return A list with `extractors` and `report`.
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
