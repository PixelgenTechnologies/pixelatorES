#' Registered Experiment Summary workflows
#'
#' Mutable package-local registry mapping workflow identifiers to zero-argument
#' functions that return extractor lists.
#'
#' @noRd
.es_data_workflow_registry <- new.env(parent = emptyenv())

#' Register an Experiment Summary workflow
#'
#' Registers a workflow identifier with a zero-argument function that returns
#' the nested extractor list used to build an `es_data` object. Extension
#' packages should call this function from their `.onLoad()` hook.
#'
#' @param name A unique workflow identifier.
#' @param extractors A zero-argument function returning a nested named list of
#'   extractor functions.
#' @param overwrite If `TRUE`, replace an existing registration for `name`.
#'   Defaults to `FALSE`.
#'
#' @return `name`, invisibly.
#'
#' @export
register_es_data_workflow <- function(name, extractors, overwrite = FALSE) {
  pixelatorR:::assert_single_value(name, "string")
  pixelatorR:::assert_class(extractors, "function")
  pixelatorR:::assert_single_value(overwrite, "bool")

  if (
    !overwrite &&
      exists(name, envir = .es_data_workflow_registry, inherits = FALSE)
  ) {
    cli_abort("Workflow {.val {name}} is already registered.")
  }

  assign(name, extractors, envir = .es_data_workflow_registry)

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

#' Get extractors for a registered workflow
#'
#' @param name A workflow identifier.
#'
#' @return A nested named list of extractor functions.
#'
#' @noRd
.get_es_data_extractors <- function(name) {
  pixelatorR:::assert_single_value(name, "string")

  if (!exists(name, envir = .es_data_workflow_registry, inherits = FALSE)) {
    cli_abort(c(
      "Unsupported {.arg workflow}: {.val {name}}.",
      "i" = "Registered workflows: {.val {list_es_data_workflows()}}."
    ))
  }

  extractor_factory <- get(
    name,
    envir = .es_data_workflow_registry,
    inherits = FALSE
  )
  extractors <- extractor_factory()
  pixelatorR:::assert_class(extractors, "list")

  return(extractors)
}
