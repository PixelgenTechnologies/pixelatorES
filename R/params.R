#' Default parameters for the ES
#'
#' @export
#'
default_params <-
  list(
    workflow = "edgelist_aggregate",
    norm_method = "CLR",
    clustering_resolution = 1,
    proximity_count_cutoff = 25,
    annotation_method = "nmf",
    mc_cores = 1,
    debug_mode = FALSE,
    test_mode = FALSE
  )

#' Print analysis parameters
#'
#' Print the analysis parameters in a table format, comparing them to the default parameters.
#'
#' @param es_data An `es_data` object containing analysis parameters.
#'
#' @return A printed table of parameters with their values and defaults.
#'
#' @export
#'
print_params <-
  function(es_data) {
    pixelatorR:::assert_class(es_data, "es_data")
    params <- es_data$params

    params[!names(params) %in% c("metadata", "sample_aliases", "sample_sheet", "data_folder")] %>%
      enframe("Parameter", "Value") %>%
      left_join(
        default_params %>%
          enframe("Parameter", "Default"),
        by = "Parameter"
      ) %>%
      left_join(
        default_params %>%
          names() %>%
          set_names() %>%
          sapply(function(x) all(params[[x]] == default_params[[x]])) %>%
          enframe("Parameter", "is_default"),
        by = "Parameter"
      ) %>%
      mutate(Default = ifelse(is_default, Default, paste("\u26A0", Default))) %>%
      select(-is_default) %>%
      style_table(caption = "Analysis parameters", interactive = FALSE, buttons = FALSE)
  }

#' Print metadata table
#'
#' Print the experiment meta data in a table format. For hashed experiments,
#' `% of pool` is included only when processed PXL data is available.
#'
#' @param es_data An `es_data` object containing the samplesheet. Processed
#'   PXL data is used when present to compute pool fractions.
#'
#' @return A printed table of sample metadata.
#'
#' @export
#'
print_metadata_table <-
  function(es_data) {
    pixelatorR:::assert_class(es_data, "es_data")
    sample_sheet <- es_data$samplesheet
    if (
      "pool" %in% names(sample_sheet) &&
        !is.null(es_data$pxl_data_processed)
    ) {
      sample_sheet <- add_pct_of_pool_to_samplesheet(
        sample_sheet,
        es_data$pxl_data_processed
      )
    }

    if (has_sample_diagnostics(es_data)) {
      flagged_aliases <- sample_diagnostic_targets(es_data)
      sample_sheet <- sample_sheet %>%
        mutate(
          Issues = ifelse(
            sample_alias %in% flagged_aliases,
            "\u26A0",
            ""
          )
        )
    }

    sample_sheet %>%
      select(
        any_of("Issues"),
        "Pool" = any_of("pool"),
        "Sample Alias" = sample_alias,
        "Sample name" = sample,
        "Condition" = condition,
        "% of pool" = any_of("pool_fraction"),
      ) %>%
      style_table(caption = "Sample settings", interactive = FALSE)
  }


#' Print Session Info table
#'
#' Print the session info in a table format.
#'
#' @return A printed table of session info.
#'
#' @export
#'
print_session_info <- function() {
  package_info() %>%
    as_tibble() %>%
    filter(attached == TRUE) %>%
    select(package, loadedversion) %>%
    style_table(caption = "Sample settings", interactive = FALSE, buttons = FALSE)
}

#' Read run information from a PXL file
#'
#' Reads the sample name and Pixelator version from the `metadata` table. The
#' panel name comes from `metadata$panel_name` in older PXL files. Newer files
#' no longer store it there, so it is built from the `panel_sources` table by
#' joining the names of all panel sources with `+`.
#'
#' @param pxl_file Path to a PXL file.
#'
#' @return A one-row tibble with `sample_name`, `version` and `panel_name`.
#'   Fields that cannot be found are `NA`.
#'
#' @noRd
.read_pxl_run_info <- function(pxl_file) {
  con <- DBI::dbConnect(duckdb::duckdb(), pxl_file, read_only = TRUE)
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)

  metadata <-
    DBI::dbGetQuery(con, "SELECT * FROM metadata")$value %>%
    RcppSimdJson::fparse()

  panel_name <- metadata$panel_name
  if (is.null(panel_name) && "panel_sources" %in% DBI::dbListTables(con)) {
    source_names <-
      DBI::dbGetQuery(con, "SELECT name FROM panel_sources ORDER BY source_id")$name
    if (length(source_names) > 0) {
      panel_name <- paste(source_names, collapse = "+")
    }
  }

  # Missing fields become NA so a single file cannot break the whole table
  as_field <- function(x) {
    if (is.null(x) || length(x) == 0) NA_character_ else as.character(x[[1]])
  }

  tibble(
    sample_name = as_field(metadata$sample_name),
    version = as_field(metadata$version),
    panel_name = as_field(panel_name)
  )
}

#' Resolve the sample alias for each PXL file in an `FSMap`
#'
#' Maps each row of `FSMap(object)` to its `sample_alias` using the cell IDs in
#' `id_map$current_id` and `object`'s cell metadata, so file-to-alias pairing
#' does not depend on the order in which PXL files were discovered or merged.
#'
#' @param object A `Seurat` object with `sample_alias` cell metadata.
#' @param fs_map A file-sample map returned by [pixelatorR::FSMap()].
#'
#' @return A character vector of sample aliases, one per row of `fs_map`.
#'
#' @noRd
.fs_map_sample_aliases <- function(object, fs_map) {
  cell_meta <- object[[]]
  if (!"sample_alias" %in% colnames(cell_meta)) {
    cli_abort(
      "{.arg es_data$pxl_data_processed} must contain {.field sample_alias} cell metadata."
    )
  }

  aliases <- vapply(
    fs_map$id_map,
    function(id_map) {
      matched_rows <- match(id_map$current_id, rownames(cell_meta))
      sample_alias <- unique(as.character(cell_meta$sample_alias[matched_rows]))
      if (length(sample_alias) != 1 || is.na(sample_alias)) {
        cli_abort(
          "Each PXL file in {.fn FSMap} must map to a single {.field sample_alias}."
        )
      }
      return(sample_alias)
    },
    character(1)
  )

  return(aliases)
}

#' Print Pixelator Version table
#'
#' Print the Pixelator version information from the pxl files.
#'
#' @param es_data An `es_data` object containing processed PXL data and the
#'   effective samplesheet.
#'
#' @return A printed table of Pixelator version information.
#' @export
#'
print_pixelator_version <-
  function(es_data) {
    pixelatorR:::assert_class(es_data, "es_data")
    object <- es_data$pxl_data_processed
    pixelatorR:::assert_class(object, "Seurat")
    sample_sheet <- es_data$effective_samplesheet

    save_fields <-
      tibble::tribble(
        ~var, ~display_name,
        "sample_alias", "Sample alias",
        "sample_name", "Sample name",
        "version", "Pixelator version",
        "panel_name", "Panel name"
      )

    fs_map <- pixelatorR::FSMap(object)
    fs_aliases <- .fs_map_sample_aliases(object, fs_map)

    sample_order <-
      if (!is.null(sample_sheet)) {
        as.character(sample_sheet$sample_alias)
      } else if (!is.null(es_data$sample_aliases)) {
        unname(es_data$sample_aliases)
      } else if (is.factor(object$sample_alias)) {
        levels(object$sample_alias)
      } else {
        fs_aliases
      }

    version_table <-
      fs_map$pxl_file %>%
      set_names(fs_aliases) %>%
      map(.read_pxl_run_info) %>%
      bind_rows(.id = "sample_alias") %>%
      arrange(match(sample_alias, sample_order)) %>%
      select(all_of(save_fields$var)) %>%
      set_names(save_fields$display_name) %>%
      style_table(
        escape = FALSE,
        search = FALSE,
        interactive = FALSE,
        buttons = FALSE
      )

    return(version_table)
  }
