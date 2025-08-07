#' Process data for clustering and dimensionality reduction
#'
#' This function processes the provided Seurat object by normalizing, scaling,
#' running PCA, and UMAP, and optionally harmonizing the data before clustering.
#'
#' @param object A Seurat object containing the data to be processed.
#' @param params A list of parameters.
#' @param annotate_cells A boolean indicating whether to annotate cells.
#'
#' @return A Seurat object with processed data, including PCA, UMAP, and clustering results.
#'
#' @export
#'
process_data <-
  function(
    object,
    params,
    annotate_cells = TRUE
  ) {
    max_dims <- ifelse(params$test_mode, 5, 20)
    npcs <- ifelse(params$test_mode, 5, 50)
    n_neighbors <- ifelse(params$test_mode, 5, 30)
    object <-
      object %>%
      NormalizeData(normalization.method = params$norm_method, margin = 2) %>%
      ScaleData() %>%
      RunPCA(npcs = npcs) %>%
      RunUMAP(reduction = "pca", dims = 1:max_dims, n.neighbors = n_neighbors)


    if (params$do_harmonize) {
      object <-
        object %>%
        RunHarmony(
          group.by.vars = params$harmonization_vars,
          plot_convergence = FALSE,
        ) %>%
        RunUMAP(
          reduction = "harmony", dims = 1:max_dims,
          reduction.name = "harmony_umap",
          n.neighbors = n_neighbors
        ) %>%
        FindNeighbors(reduction = "harmony", dims = 1:max_dims) %>%
        FindClusters(random.seed = 1, resolution = params$clustering_resolution)
    } else {
      object <-
        object %>%
        FindNeighbors(reduction = "pca", dims = 1:max_dims) %>%
        FindClusters(random.seed = 1, resolution = params$clustering_resolution)
    }

    if (annotate_cells) {
      reference <- read_annotation_reference()
      reference$l1_annotation <- reference$celltype.l1
      reference$l2_annotation <- reference$celltype.l2
      object <-
        object %>%
        AnnotateCells(
          reference = reference,
          reference_groups = c("l1_annotation", "l2_annotation"),
          reference_assay = "PNA",
          query_assay = "PNA",
          summarize_by_column = "seurat_clusters",
          normalization_method = ifelse(params$annotation_method == "Seurat", "LogNormalize", "CLR"),
          method = params$annotation_method
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
#' @param sample_levels Optional vector of sample levels to order the samples in the plots.
#'
#' @return A data frame of filtered proximity scores with additional metadata.
#'
#' @export
#'
filter_proximity_scores <- function(
  pg_data_processed,
  params,
  sample_levels = NULL
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
      ))
    )

  proximity_scores <- set_sample_levels(proximity_scores, sample_levels)

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
