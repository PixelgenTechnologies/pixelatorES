#' Default parameters for the ES
#'
#' @export
#'
default_params <-
  list(
    control_markers = c("mIgG1", "mIgG2a", "mIgG2b"),
    do_harmonize = FALSE,
    harmonization_vars = "condition",
    norm_method = "CLR",
    clustering_resolution = 1,
    proximity_count_cutoff = 25,
    annotation_method = "Seurat",
    mc_cores = 1,
    debug_mode = FALSE,
    test_mode = FALSE,
  )

#' Print analysis parameters
#'
#' Print the analysis parameters in a table format, comparing them to the default parameters.
#'
#' @param params A list of parameters used in the analysis.
#'
#' @return A printed table of parameters with their values and defaults.
#'
#' @export
#'
print_params <-
  function(params) {
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
#' @param sample_sheet A sample sheet.
#'
#' @return A printed table of sample metadata.
#'
#' @export
#'
print_metadata_table <-
  function(sample_sheet) {
    sample_sheet %>%
      select(
        "Sample Alias" = sample_alias,
        "Sample name" = sample,
        "Condition" = condition
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
