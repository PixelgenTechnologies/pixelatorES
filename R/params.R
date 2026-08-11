#' Default parameters for the ES
#'
#' @export
#'
default_params <-
  list(
    workflow = "amplicon_demux",
    control_markers = c("mIgG1", "mIgG2a", "mIgG2b"),
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
#' Print the experiment meta data in a table format.
#'
#' @param es_data An `es_data` object containing the samplesheet and processed
#'   PXL data.
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

    fs_map$pxl_file %>%
      set_names(sample_sheet$sample_alias) %>%
      sapply(function(x) {
        # Open connection
        con <-
          DBI::dbConnect(
            duckdb::duckdb(),
            x,
            read_only = TRUE
          )

        # Get metadata
        metadata <-
          DBI::dbGetQuery(
            con,
            "SELECT * FROM metadata"
          )

        # Close connection
        DBI::dbDisconnect(con, shutdown = TRUE)

        metadata_parsed <-
          metadata$value %>%
          RcppSimdJson::fparse()

        return(metadata_parsed)
      }) %>%
      t() %>%
      as_tibble(rownames = "sample_alias") %>%
      select(all_of(save_fields$var)) %>%
      unnest(cols = everything()) %>%
      set_names(save_fields$display_name) %>%
      style_table(
        escape = FALSE,
        search = FALSE,
        interactive = FALSE,
        buttons = FALSE
      )
  }
