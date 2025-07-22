#' Annotate cells in a Seurat object using a reference dataset
#'
#' This function annotates cells in a Seurat object using a reference dataset.
#'
#' @param object A Seurat object containing the query data.
#' @param reference A Seurat object containing the reference data.
#' @param summarize_by_column Optional. A column name in the object metadata to summarize annotations by.
#' @param reference_assay The assay in the reference object to use for annotation. Default is "PNA".
#' @param query_assay The assay in the query object to use for annotation. Default is "PNA".
#' @param normalization_method The normalization method to use. Default is "CLR".
#' @param reduction The reduction method to use for finding transfer anchors. Default is "cca".
#' @param ... Additional arguments passed to `FindTransferAnchors` and `TransferData`.
#'
#' @return A Seurat object with cell annotations added to the metadata.
#'
#' @export
#'
annotate_cells <-
  function(
             object,
             reference,
             summarize_by_column = NULL,
             reference_assay = "PNA",
             query_assay = "PNA",
             normalization_method = "LogNormalize",
             reduction = "cca",
             ...) {
    pixelatorR:::assert_class(object, "Seurat")
    pixelatorR:::assert_class(reference, "Seurat")
    pixelatorR:::assert_class(summarize_by_column, "character",
      allow_null = TRUE
    )
    pixelatorR:::assert_class(reference_assay, "character")
    pixelatorR:::assert_class(query_assay, "character")
    pixelatorR:::assert_class(normalization_method, "character")
    normalization_method <-
      match.arg(
        normalization_method,
        c("LogNormalize", "SCT")
      )
    pixelatorR:::assert_class(reduction, "character")

    stopifnot(
      `'object' must be a Seurat object.` = inherits(
        object,
        "Seurat"
      ),
      `'reference' must be a Seurat object.` = inherits(
        reference,
        "Seurat"
      ),
      `'summarize_by_column' must be NULL or present in the object metadata.` = is.null(summarize_by_column) ||
        summarize_by_column %in% colnames(object[[]])
    )

    tmp <-
      object %>%
      NormalizeData(
        assay = query_assay,
        normalization.method = normalization_method
      )

    reference <-
      reference %>%
      NormalizeData(
        assay = reference_assay,
        normalization.method = normalization_method
      )

    feats <-
      intersect(rownames(tmp), rownames(reference))

    anchors <- Seurat::FindTransferAnchors(
      reference = reference,
      query = tmp,
      reference.assay = reference_assay,
      query.assay = query_assay,
      normalization.method = normalization_method,
      reduction = reduction,
      features = feats,
      l2.norm = TRUE
    )

    predictions <- Seurat::TransferData(
      anchorset = anchors,
      refdata = list(
        celltype.l1 = reference@meta.data$celltype.l1,
        celltype.l2 = reference@meta.data$celltype.l2
      ),
      weight.reduction = reduction
    )

    object <- SeuratObject::AddMetaData(
      object = object,
      metadata = data.frame(
        l2_annotation = predictions$celltype.l2[, "predicted.id"],
        l1_annotation = predictions$celltype.l1[, "predicted.id"]
      )
    )


    if (!is.null(summarize_by_column)) {
      get_annotation_summary <-
        . %>%
        group_by(cluster, anno) %>%
        count() %>%
        group_by(cluster) %>%
        top_n(1, n) %>%
        slice_head(n = 1) %>%
        select(cluster, anno) %>%
        deframe()


      annotation_summary <-
        FetchData(object, c(summarize_by_column, "l1_annotation", "l2_annotation")) %>%
        rename(cluster = !!sym(summarize_by_column)) %>%
        as_tibble(rownames = "id") %>%
        mutate(
          l1_annotation_summary = get_annotation_summary(rename(.,
            anno = l1_annotation
          ))[cluster],
          l2_annotation_summary = get_annotation_summary(rename(.,
            anno = l2_annotation
          ))[cluster]
        ) %>%
        select(id, l1_annotation_summary, l2_annotation_summary) %>%
        column_to_rownames("id")

      object <- SeuratObject::AddMetaData(
        object = object,
        metadata = annotation_summary
      )
    }
    return(object)
  }


#' Read cell annotation reference
#'
#' This function reads the cell annotation reference data from the package's extdata directory.
#'
#' @return A Seurat object containing the cell annotation and counts data.
#'
#' @export
#'
read_annotation_reference <-
  function() {
    cell_annotation <-
      system.file("extdata", "PBMC_cell_annotation.csv",
        package = "pixelatorESinternal"
      ) %>%
      read_csv(col_types = "c") %>%
      column_to_rownames("sample_component")

    counts <-
      system.file("extdata", "PBMC_cell_counts.csv",
        package = "pixelatorESinternal"
      ) %>%
      read_csv(col_types = cols(
        sample_component = col_character(),
        .default = col_integer()
      )) %>%
      column_to_rownames("sample_component")

    counts <-
      counts[rownames(cell_annotation), ] %>%
      t() %>%
      as("dgCMatrix")

    # Create Seurat object

    CreateSeuratObject(
      counts = counts,
      assay = "PNA",
      meta.data = cell_annotation
    )
  }
