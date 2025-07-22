#' Process data for clustering and dimensionality reduction
#'
#' This function processes the provided Seurat object by normalizing, scaling,
#' running PCA, and UMAP, and optionally harmonizing the data before clustering.
#'
#' @param object A Seurat object containing the data to be processed.
#' @param norm_method A string specifying the normalization method to use (default is "CLR").
#' @param clustering_resolution A numeric value specifying the resolution for clustering (default is 1).
#' @param do_harmonize A boolean indicating whether to perform data harmonization (default is FALSE).
#' @param harmonization_vars A character vector specifying the variables to use for harmonization (default is NULL).
#' @param annotate_cells A boolean indicating whether to annotate cells using a reference dataset (default is TRUE).
#' @param test_mode A boolean indicating whether to run in test mode (default is FALSE).
#'
#' @return A Seurat object with processed data, including PCA, UMAP, and clustering results.
#'
#' @export
#'
process_data <-
  function(
    object,
    norm_method = "CLR",
    clustering_resolution = 1,
    do_harmonize = FALSE,
    harmonization_vars = NULL,
    annotate_cells = TRUE,
    test_mode = FALSE
  ) {
    max_dims <- ifelse(test_mode, 5, 20)
    npcs <- ifelse(test_mode, 5, 50)
    n_neighbors <- ifelse(test_mode, 5, 30)
    object <-
      object %>%
      NormalizeData(normalization.method = norm_method, margin = 2) %>%
      ScaleData() %>%
      RunPCA(npcs = npcs) %>%
      RunUMAP(reduction = "pca", dims = 1:max_dims, n.neighbors = n_neighbors)


    if (do_harmonize) {
      object <-
        object %>%
        RunHarmony(
          group.by.vars = harmonization_vars,
          plot_convergence = FALSE,
        ) %>%
        RunUMAP(
          reduction = "harmony", dims = 1:max_dims,
          reduction.name = "harmony_umap",
          n.neighbors = n_neighbors
        ) %>%
        FindNeighbors(reduction = "harmony", dims = 1:max_dims) %>%
        FindClusters(random.seed = 1, resolution = clustering_resolution)
    } else {
      object <-
        object %>%
        FindNeighbors(reduction = "pca", dims = 1:max_dims) %>%
        FindClusters(random.seed = 1, resolution = clustering_resolution)
    }

    if (annotate_cells && !test_mode) {
      object <-
        object %>%
        pixelatorES::annotate_cells(
          reference = read_annotation_reference(),
          reference_assay = "PNA",
          query_assay = "PNA",
          summarize_by_column = "seurat_clusters"
        )
    }

    return(object)
  }


#' Filter Proximity Scores
#'
#' This function filters the proximity scores from the `pg_data_processed` object.
#'
#' @param pg_data_processed A Seurat object.
#' @param params A list of parameters for filtering, including control markers.
#' @param metadata A data frame containing metadata for the samples.
#'
#' @return A data frame of filtered proximity scores with additional metadata.
#'
#' @export
#'
filter_proximity_scores <- function(
  pg_data_processed,
  params,
  metadata
) {
  # Set proximity score cutoff
  background_threshold_pct <-
    median(pg_data_processed$isotype_fraction) +
    2 * mad(pg_data_processed$isotype_fraction)

  # Fall back to a small cutoff if the background threshold is too low
  background_threshold_pct <- max(background_threshold_pct, 0.001)

  proximity_scores <-
    ProximityScores(pg_data_processed,
      add_marker_counts = TRUE, add_marker_proportions = TRUE,
      meta_data_columns = c(
        "sample_alias", "condition",
        "seurat_clusters", "l1_annotation_summary",
        "condition"
      ),
      lazy = TRUE
    ) %>%
    rename(sample_component = component) %>%
    FilterProximityScores(
      background_threshold_pct = background_threshold_pct
    ) %>%
    collect() %>%
    mutate(
      marker_1 = factor(marker_1, order_cd_markers(
        unique(marker_1),
        params$control_markers
      )),
      marker_2 = factor(marker_2, order_cd_markers(
        unique(marker_2),
        params$control_markers
      )),
      sample_alias = factor(sample_alias, levels = metadata$sample_alias),
      condition = factor(condition, levels = unique(metadata$condition))
    )

  return(proximity_scores)
}


#' Process Proximity Scores
#'
#' This function processes the proximity scores to prepare them for plotting.
#'
#' @param pg_data_processed A Seurat object containing processed data.
#' @param proximity_scores A data frame of proximity scores.
#'
#' @return A data frame.
#'
#' @export
#'
process_proximity_scores <- function(
  pg_data_processed,
  proximity_scores
) {
  plot_contrasts <-
    c(
      "B2M" = "HLA-ABC",

      # Sgnalling microclusters
      "CD45" = "CD45RA",
      "CD45" = "CD45RB",
      "CD45" = "CD45RO",
      "CD11a" = "CD18",
      "CD29" = "CD49D",
      "CD11b" = "CD18",
      "CD11c" = "CD18",

      # Isoforms
      "HLA-DR" = "HLA-DR-DP-DQ",

      # T cell
      "CD3e" = "TCRab",
      "CD2" = "CD58",

      # B cell
      "CD19" = "CD21",
      "CD19" = "CD81",
      "CD79a" = "IgM",
      "CD79a" = "IgD",

      # NK cell
      "CD159a" = "CD94",

      # Tetraspanins
      "CD81" = "CD82",
      "CD53" = "CD82",

      # Complement cascade
      "CD55" = "CD59",
      "CD21" = "CD35"
    ) %>%
    enframe("marker_1", "marker_2") %>%
    unite("contrast", marker_1, marker_2, sep = "/", remove = FALSE)

  processed_data <-
    proximity_scores %>%
    filter(l1_annotation_summary %in% c(
      "CD4 T", "CD8 T",
      "NK", "Mono", "B"
    )) %>%
    inner_join(plot_contrasts) %>%
    left_join(
      FetchData(pg_data_processed, c(
        "sample_alias",
        "l1_annotation_summary"
      )) %>%
        group_by_all() %>%
        count(name = "n_cells") %>%
        group_by(sample_alias) %>%
        mutate(n_cells_tot = sum(n_cells)),
      by = c("sample_alias", "l1_annotation_summary")
    ) %>%
    group_by(contrast)

  return(processed_data)
}


#' Summarize proximity scores per sample and condition
#'
#' This function summarizes the proximity per sample and condition.
#'
#' @param proximity_scores A data frame of proximity scores.
#' @param test_mode A logical indicating whether to run in test mode (default is FALSE).
#'
#' @return A data frame summarizing the proximity scores per sample and condition.
#'
#' @export
#'
summarize_colocalization_scores_per_sample <- function(
  proximity_scores,
  test_mode = FALSE
) {
  summarized_data <-
    proximity_scores %>%
    rename(component = sample_component) %>%
    SummarizeProximityScores(
      proximity_metric = "join_count_z",
      detailed = TRUE,
      group_vars = c("sample_alias", "condition")
    ) %>%
    rowwise() %>%
    mutate(mean_log2_ratio = mean(log2(
      pmax(unlist(join_count_list), 1) /
        pmax(unlist(join_count_expected_mean_list), 1)
    ))) %>%
    {
      if (test_mode) {
        .
      } else {
        filter(., n_cells_detected >= 10)
      }
    } %>%
    ungroup() %>%
    select(
      sample_alias, condition, marker_1, marker_2,
      mean_join_count_z, mean_log2_ratio
    ) %>%
    group_by(sample_alias)

  return(summarized_data)
}


#' Summarize proximity scores per celltype and condition
#'
#' This function summarizes the proximity per celltype and condition.
#'
#' @param proximity_scores A data frame of proximity scores.
#' @param test_mode A logical indicating whether to run in test mode (default is FALSE).
#'
#' @return A data frame summarizing the proximity scores per celltype and condition.
#'
#' @export
#'
summarize_colocalization_scores_per_celltype <- function(
  proximity_scores,
  test_mode = FALSE
) {
  summarized_data <-
    proximity_scores %>%
    rename(component = sample_component) %>%
    SummarizeProximityScores(
      proximity_metric = "join_count_z",
      detailed = TRUE,
      group_vars = c(
        "condition",
        "l1_annotation_summary"
      )
    ) %>%
    rowwise() %>%
    mutate(mean_log2_ratio = mean(log2(
      pmax(unlist(join_count_list), 1) /
        pmax(unlist(join_count_expected_mean_list), 1)
    ))) %>%
    {
      if (test_mode) {
        .
      } else {
        filter(., n_cells_detected >= 10)
      }
    } %>%
    ungroup() %>%
    select(
      condition, l1_annotation_summary, marker_1, marker_2,
      mean_join_count_z, mean_log2_ratio
    ) %>%
    group_by(celltype = l1_annotation_summary, condition)

  return(summarized_data)
}
