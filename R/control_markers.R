#' Get control markers from the panel metadata
#'
#' Detects which antibodies are control antibodies (e.g. isotype controls) from
#' the `control` column in the feature metadata of the default assay. This
#' column is written by Pixelator from the panel definition, so the control
#' markers follow the panel and do not need to be specified by the user. The
#' column can be logical or contain `"yes"`/`"no"`.
#'
#' @param object A Seurat object with a `PNAAssay` or `PNAAssay5` default assay,
#'   for example from `pixelatorR::ReadPNA_Seurat()`.
#'
#' @return A character vector of control marker names. The vector is empty if
#'   the panel has no control markers.
#'
#' @export
#'
get_control_markers <- function(object) {
  pixelatorR:::assert_class(object, "Seurat")

  assay <- object[[SeuratObject::DefaultAssay(object)]]
  feature_meta <-
    if (methods::is(assay, "Assay5")) {
      assay[[]]
    } else {
      methods::slot(assay, "meta.features")
    }

  if (!"control" %in% colnames(feature_meta)) {
    cli_abort(
      c(
        "Cannot detect control markers: the feature metadata has no {.field control} column.",
        "i" = "The {.field control} column is written by Pixelator from the panel definition.",
        "i" = "Use pixelator and pixelatorR versions that include it in the PXL file."
      )
    )
  }

  is_control <-
    if (is.logical(feature_meta$control)) {
      feature_meta$control %in% TRUE
    } else {
      tolower(as.character(feature_meta$control)) %in% c("yes", "true")
    }
  control_markers <- rownames(feature_meta)[is_control]

  return(control_markers)
}

#' Detect control markers across samples
#'
#' Takes the union of the control markers of each sample, so that samples with
#' different panels are all covered.
#'
#' @param pxl_data A list of Seurat objects, one per sample.
#'
#' @return A sorted, unique character vector of control marker names.
#'
#' @noRd
.get_control_markers_union <- function(pxl_data) {
  control_markers <-
    pxl_data %>%
    map(get_control_markers) %>%
    unlist(use.names = FALSE) %>%
    as.character() %>%
    unique() %>%
    sort()

  return(control_markers)
}
