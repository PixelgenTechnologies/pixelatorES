#' Technologies accepted by the Experiment Summary
#'
#' These are the same values as the pipeline `--technology` parameter.
#' `proxiome-v1` and `nonhashed_samples` are unhashed: a samplesheet `pool`
#' column is ignored. `proxiome-v2` and `hashed_samples` are hashed: a
#' samplesheet `pool` column is required.
#'
#' @noRd
.es_technology_unhashed <- c("proxiome-v1", "nonhashed_samples")

#' @noRd
.es_technology_hashed <- c("proxiome-v2", "hashed_samples")

#' @noRd
.es_technologies <- c(.es_technology_unhashed, .es_technology_hashed)

#' Check a technology value
#'
#' The Experiment Summary does not infer hashing from the samplesheet or from
#' the files in the data folder. The caller sets `technology` up front, using
#' one of the pipeline values.
#'
#' @param technology Technology name supplied with the report parameters.
#'
#' @return `technology` when it is one of `proxiome-v1`, `nonhashed_samples`,
#'   `proxiome-v2`, or `hashed_samples`.
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
        "x" = "{.arg technology} must be one of {.val { .es_technologies }}.",
        "i" = paste(
          "{.val { .es_technology_unhashed }} ignore a samplesheet",
          "{.field pool} column.",
          "{.val { .es_technology_hashed }} require one."
        )
      )
    )
  }

  return(as.character(technology))
}

#' Whether a technology uses sample hashing
#'
#' @param technology Technology name. See [validate_es_technology()].
#'
#' @return `TRUE` for `proxiome-v2` and `hashed_samples`. `FALSE` for
#'   `proxiome-v1` and `nonhashed_samples`.
#'
#' @export
#'
is_hashed_technology <- function(technology) {
  return(validate_es_technology(technology) %in% .es_technology_hashed)
}

#' Apply a technology to a samplesheet
#'
#' @param sample_sheet A samplesheet, as returned by [read_samplesheet()].
#' @param technology A pipeline technology name. See [validate_es_technology()].
#'
#' @return The samplesheet used for the rest of the run. For an unhashed
#'   technology the `pool` column is removed when present. For a hashed
#'   technology the `pool` column is kept and must contain a non-empty value
#'   on every row.
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
          "Technology {.val {technology}} requires a",
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
          "Technology {.val {technology}} requires a non-empty",
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
