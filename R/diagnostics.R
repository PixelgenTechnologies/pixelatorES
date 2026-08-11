#' Convert Experiment Summary diagnostics to a table
#'
#' Converts the diagnostics recorded while building an `es_data` object to a
#' tibble with one row per diagnostic.
#'
#' @param es_data An `es_data` object.
#'
#' @return A tibble with the columns `type`, `target`, and `message`.
#'
#' @export
#'
diagnostics_to_tibble <- function(es_data) {
  pixelatorR:::assert_class(es_data, "es_data")

  if (length(es_data$diagnostics) == 0) {
    return(tibble(
      type = character(),
      target = character(),
      message = character()
    ))
  }

  diagnostics <- bind_rows(es_data$diagnostics) %>%
    select(type, target, message)

  return(diagnostics)
}

#' Label a diagnostic type
#'
#' Translates the internal diagnostic types to labels shown in the report.
#'
#' @param type A character vector of diagnostic types.
#'
#' @return A character vector of labels.
#'
#' @noRd
#'
.diagnostic_type_label <- function(type) {
  labels <- c(
    pxl_load = "PXL loading",
    qc_load = "QC loading",
    extractor = "Analysis step"
  )

  return(unname(labels[type]))
}

#' Find samples affected by diagnostics
#'
#' Finds sample aliases targeted by PXL or QC loading diagnostics. Diagnostics
#' targeting a pool affect every sample assigned to that pool.
#'
#' @param es_data An `es_data` object.
#'
#' @return A character vector of affected sample aliases.
#'
#' @noRd
#'
sample_diagnostic_targets <- function(es_data) {
  pixelatorR:::assert_class(es_data, "es_data")

  sample_sheet <- es_data$samplesheet
  if (is.null(sample_sheet) || nrow(sample_sheet) == 0) {
    return(character())
  }

  diagnostics <- diagnostics_to_tibble(es_data) %>%
    filter(type %in% c("pxl_load", "qc_load"))

  if (nrow(diagnostics) == 0) {
    return(character())
  }

  sample_aliases <- as.character(sample_sheet$sample_alias)
  affected_aliases <- diagnostics$target[
    diagnostics$target %in% sample_aliases
  ]

  if ("pool" %in% names(sample_sheet)) {
    affected_pools <- diagnostics$target[
      diagnostics$target %in% as.character(sample_sheet$pool)
    ]
    pool_aliases <- sample_sheet %>%
      filter(as.character(pool) %in% affected_pools) %>%
      pull(sample_alias) %>%
      as.character()
    affected_aliases <- c(affected_aliases, pool_aliases)
  }

  return(unique(affected_aliases))
}

#' Check whether samples have diagnostics
#'
#' Checks whether any PXL or QC loading diagnostic targets a sample or its
#' pool.
#'
#' @param es_data An `es_data` object.
#'
#' @return `TRUE` when at least one sample is affected; otherwise `FALSE`.
#'
#' @export
#'
has_sample_diagnostics <- function(es_data) {
  return(length(sample_diagnostic_targets(es_data)) > 0)
}

#' Format sample diagnostics for an Experiment Summary
#'
#' Formats sample- and pool-targeted loading diagnostics as Markdown lines for
#' display on the Samples page.
#'
#' @param es_data An `es_data` object.
#'
#' @return A single Markdown string, or `NULL` when no samples are affected.
#'
#' @export
#'
format_sample_diagnostics_summary <- function(es_data) {
  pixelatorR:::assert_class(es_data, "es_data")

  sample_sheet <- es_data$samplesheet
  if (is.null(sample_sheet) || nrow(sample_sheet) == 0) {
    return(NULL)
  }

  diagnostics <- diagnostics_to_tibble(es_data) %>%
    filter(type %in% c("pxl_load", "qc_load"))

  if ("pool" %in% names(sample_sheet)) {
    valid_targets <- c(
      as.character(sample_sheet$sample_alias),
      as.character(sample_sheet$pool)
    )
  } else {
    valid_targets <- as.character(sample_sheet$sample_alias)
  }

  diagnostics <- diagnostics %>%
    filter(target %in% valid_targets)

  if (nrow(diagnostics) == 0) {
    return(NULL)
  }

  lines <- paste0(
    "- **",
    diagnostics$target,
    "** (",
    .diagnostic_type_label(diagnostics$type),
    "): ",
    diagnostics$message
  )

  return(paste(lines, collapse = "\n"))
}

#' Format a sample diagnostics callout for an Experiment Summary
#'
#' Formats the sample- and pool-targeted loading diagnostics as a red Quarto
#' callout listing every affected sample, so that a partially loaded experiment
#' is immediately visible on the Samples page.
#'
#' @param es_data An `es_data` object.
#'
#' @return A single Markdown string holding the callout, or `NULL` when no
#'   samples are affected.
#'
#' @export
#'
format_sample_diagnostics_callout <- function(es_data) {
  summary_lines <- format_sample_diagnostics_summary(es_data)

  if (is.null(summary_lines)) {
    return(NULL)
  }

  callout <- paste0(
    '::: {.callout-important title="Sample loading issues"}\n',
    "Some samples could not be fully loaded, and the metrics in this report ",
    "are therefore incomplete. Affected samples are marked with a warning ",
    "symbol in the table below.\n\n",
    summary_lines,
    "\n\nSee the Diagnostics section under Run info for the complete list.\n",
    ":::\n"
  )

  return(callout)
}

#' Print Experiment Summary diagnostics
#'
#' Prints all diagnostics recorded while building an `es_data` object as a
#' styled table.
#'
#' @param es_data An `es_data` object.
#'
#' @return A `datatables` HTML widget containing all diagnostics.
#'
#' @export
#'
print_diagnostics_table <- function(es_data) {
  diagnostics <- diagnostics_to_tibble(es_data) %>%
    mutate(type = .diagnostic_type_label(type)) %>%
    rename(
      "Type" = type,
      "Target" = target,
      "Message" = message
    ) %>%
    style_table(interactive = FALSE, escape = TRUE)

  return(diagnostics)
}
