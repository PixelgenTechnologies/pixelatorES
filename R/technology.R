#' Proxiome technologies accepted by the Experiment Summary
#'
#' `proxiome-v1` is an unhashed run: a samplesheet `pool` column is ignored.
#' `proxiome-v2` is a hashed run: a samplesheet `pool` column is required.
#'
#' @noRd
.es_technology_unhashed <- "proxiome-v1"

#' @noRd
.es_technology_hashed <- "proxiome-v2"

#' @noRd
.es_technologies <- c(.es_technology_unhashed, .es_technology_hashed)

#' Check a Proxiome technology value
#'
#' The Experiment Summary does not infer hashing from the samplesheet or from
#' the files in the data folder. The caller sets `technology` up front.
#'
#' @param technology Technology name supplied with the report parameters.
#'
#' @return `technology` when it is `proxiome-v1` or `proxiome-v2`.
#'
#' @export
#'
validate_es_technology <- function(technology) {
  technology_ok <-
    !is.null(technology) &&
      length(technology) == 1 &&
      !is.na(technology) &&
      as.character(technology) %in% .es_technologies

  if (!technology_ok) {
    cli_abort(
      c(
        "x" = paste(
          "{.arg technology} must be {.val { .es_technology_unhashed }}",
          "or {.val { .es_technology_hashed }}."
        ),
        "i" = paste(
          "{.val { .es_technology_unhashed }} ignores a samplesheet {.field pool} column.",
          "{.val { .es_technology_hashed }} requires one."
        )
      )
    )
  }

  return(as.character(technology))
}

#' Whether a technology uses sample hashing
#'
#' @param technology Technology name, `proxiome-v1` or `proxiome-v2`.
#'
#' @return `TRUE` for `proxiome-v2`, `FALSE` for `proxiome-v1`.
#'
#' @export
#'
is_hashed_technology <- function(technology) {
  return(identical(
    validate_es_technology(technology),
    .es_technology_hashed
  ))
}

#' Apply a technology to a samplesheet
#'
#' @param sample_sheet A samplesheet, as returned by [read_samplesheet()].
#' @param technology `proxiome-v1` or `proxiome-v2`. See
#'   [validate_es_technology()].
#'
#' @return The samplesheet used for the rest of the run. For `proxiome-v1` the
#'   `pool` column is removed when present. For `proxiome-v2` the `pool`
#'   column is kept and must contain a non-empty value on every row.
#'
#' @export
#'
apply_es_technology <- function(sample_sheet, technology) {
  pixelatorR:::assert_class(sample_sheet, "tbl_df")
  technology <- validate_es_technology(technology)

  if (!is_hashed_technology(technology)) {
    sample_sheet$pool <- NULL
    return(sample_sheet)
  }

  if (!"pool" %in% names(sample_sheet)) {
    cli_abort(
      c(
        "x" = paste(
          "Technology {.val { .es_technology_hashed }} requires a",
          "{.field pool} column in the samplesheet."
        )
      )
    )
  }

  pool_values <- as.character(sample_sheet$pool)
  pool_missing <- is.na(pool_values) | !nzchar(pool_values)
  if (any(pool_missing)) {
    cli_abort(
      c(
        "x" = paste(
          "Technology {.val { .es_technology_hashed }} requires a non-empty",
          "{.field pool} value for every sample."
        )
      )
    )
  }

  # Pool file names and sample file names are matched on the same basename,
  # so a pool id that is also a sample id would hide that sample's files.
  pool_overlaps_sample <-
    c(sample_sheet$sample, sample_sheet$sample_alias) %in% pool_values
  if (any(pool_overlaps_sample)) {
    cli_abort(
      c(
        "x" = paste(
          "The {.field pool} column cannot contain values that are also",
          "present in {.field sample} or {.field sample_alias}."
        )
      )
    )
  }

  return(sample_sheet)
}
