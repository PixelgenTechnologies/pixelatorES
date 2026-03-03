#' Cell types to display
#'
#' This vector contains the cell types to be displayed in the ES, when filtering is applied.
#'
#' @export
#'
displayed_cell_types <-
  c("T", "NK", "B", "Mono & DC", "Platelets")

#' Cell type merging
#'
#' This vector contains the mapping of original cell types to merged cell types, which is used by the ES to merge cell
#' types for display purposes. This allows for a more concise and interpretable visualization.
#'
#' @export
#'
merged_cell_types <-
  c(
    "CD4 T" = "T",
    "CD8 T" = "T",
    "Other T" = "T",
    "NK" = "NK",
    "Mono" = "Mono & DC",
    "B" = "B",
    "Platelets" = "Platelets"
  )

#' Merge cell types
#'
#' This function takes a Seurat object and merges the cell types based on the provided mapping. It adds a new column to
#' the metadata with the merged cell types.
#'
#' @param object A Seurat object containing the cell type annotations in the metadata.
#' @param celltype_col The name of the column in the metadata that contains the original cell type annotations.
#' Default is "l1_annotation_summary".
#' @param merged_col The name of the new column to be added to the metadata that will contain the merged cell type
#' annotations. Default is "celltype".
#' @param mapping A named vector that maps original cell types to merged cell types. Default is the `merged_cell_types`
#' vector.
#'
#' @return A Seurat object with the merged cell types added to the metadata.
#'
#' @export
#'
merge_cell_types <-
  function(
    object,
    celltype_col = "l1_annotation_summary",
    merged_col = "celltype",
    mapping = merged_cell_types
  ) {
    object[[]] <-
      object[[]] %>%
      as_tibble(rownames = "cell_id") %>%
      mutate(
        !!sym(merged_col) := recode(!!sym(celltype_col), !!!mapping)
      ) %>%
      column_to_rownames("cell_id")

    return(object)
  }
