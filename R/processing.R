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
    stats::quantile(pg_data_processed$isotype_fraction, 0.9)

  # Fall back to a small cutoff if the background threshold is too low
  background_threshold_pct <- max(background_threshold_pct, 0.0025)

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

#' Complete Proximity Scores
#'
#' This function completes the proximity scores by ensuring that each marker pair has a score for each sample component,
#' filling missing values with 0.
#'
#' @param proximity_scores A data frame of proximity scores.
#' @param only_self A boolean indicating whether to filter for self-comparisons only (default is TRUE).
#'
#' @return A data frame of completed proximity scores with missing values filled.
#'
#' @export
#'
complete_proximity_scores <-
  function(
    proximity_scores,
    only_self = TRUE
  ) {
    pixelatorR:::assert_class(proximity_scores, "tbl_df")
    pixelatorR:::assert_single_value(only_self, "bool")

    if (only_self) {
      proximity_scores <-
        proximity_scores %>%
        filter(marker_1 == marker_2)
    }

    proximity_scores %>%
      complete(
        nesting(sample_component, sample_alias, condition, seurat_clusters, l1_annotation_summary),
        nesting(marker_1, marker_2),
        fill = list(log2_ratio = 0)
      )
  }

#' Summarize proximity scores per sample and condition
#'
#' This function summarizes the proximity per sample and condition.
#'
#' @param proximity_scores A data frame of proximity scores.
#' @param plot_markers A character vector of markers to plot (default is NULL, which means all markers will be used).
#' @param test_mode A logical indicating whether to run in test mode (default is FALSE).
#'
#' @return A data frame summarizing the proximity scores per sample and condition.
#'
#' @export
#'
summarize_colocalization_scores_per_sample <- function(
  proximity_scores,
  plot_markers = NULL,
  test_mode = FALSE
) {
  pixelatorR:::assert_class(proximity_scores, "tbl_df")
  pixelatorR:::assert_single_value(test_mode, "bool")
  pixelatorR:::assert_vector(plot_markers, "character", allow_null = TRUE)

  if (!is.null(plot_markers)) {
    proximity_scores <-
      proximity_scores %>%
      filter(marker_1 %in% plot_markers & marker_2 %in% plot_markers)
  }

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
      mean_join_count_z, mean_log2_ratio, pct_detected
    ) %>%
    complete(
      nesting(sample_alias, condition),
      nesting(marker_1, marker_2),
      fill = list(mean_join_count_z = 0, mean_log2_ratio = 0, pct_detected = 0)
    )

  return(summarized_data)
}


#' Summarize proximity scores per celltype and condition
#'
#' This function summarizes the proximity per celltype and condition.
#'
#' @param proximity_scores A data frame of proximity scores.
#' @param plot_markers A character vector of markers to plot (default is NULL, which means all markers will be used).
#' @param n_markers An integer specifying the number of top markers to select (default is 40).
#' @param test_mode A logical indicating whether to run in test mode (default is FALSE).
#'
#' @return A data frame summarizing the proximity scores per celltype and condition.
#'
#' @export
#'
summarize_colocalization_scores_per_celltype <- function(
  proximity_scores,
  plot_markers = NULL,
  n_markers = 40,
  test_mode = FALSE
) {
  pixelatorR:::assert_class(proximity_scores, "tbl_df")
  pixelatorR:::assert_single_value(test_mode, "bool")
  pixelatorR:::assert_single_value(n_markers, "numeric")
  pixelatorR:::assert_vector(plot_markers, "character", allow_null = TRUE)


  if (!is.null(plot_markers)) {
    proximity_scores <-
      proximity_scores %>%
      filter(marker_1 %in% plot_markers & marker_2 %in% plot_markers)
  }

  summarized_data <-
    proximity_scores %>%
    rename(component = sample_component) %>%
    SummarizeProximityScores(
      proximity_metric = "join_count_z",
      detailed = TRUE,
      group_vars = c(
        "sample_alias", "condition",
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
      sample_alias, condition, l1_annotation_summary, marker_1, marker_2,
      mean_join_count_z, mean_log2_ratio, pct_detected
    ) %>%
    complete(
      nesting(sample_alias, condition, l1_annotation_summary),
      nesting(marker_1, marker_2),
      fill = list(mean_join_count_z = 0, mean_log2_ratio = 0, pct_detected = 0)
    )

  return(summarized_data)
}

#' Run ANOVA for each marker
#'
#' This function runs an ANOVA for each marker in the dataset.
#'
#' @param data A data frame containing columns `x`, `sample`, and `cluster`.
#' @param vars A character vector of length 1 or 2 specifying the variables to use in the ANOVA.
#' @param response A string specifying the response variable (default is "x").
#' @param formula_str A formula specifying the model to use in the ANOVA (optional).
#' @param tidy Logical indicating whether to return a tidy data frame with ANOVA results.
#'
#' @return A model or data frame with ANOVA results (if `tidy = TRUE`).
#'
#' @noRd
#'
.run_anova <- function(
  data,
  vars = NULL,
  response = "x",
  formula_str = NULL,
  tidy = TRUE
) {
  pixelatorR:::assert_class(data, "data.frame")
  pixelatorR:::assert_vector(vars, "character", n = 1, allow_null = TRUE)
  pixelatorR:::assert_single_value(formula_str, "string", allow_null = TRUE)
  pixelatorR:::assert_max_length(vars, 2, allow_null = TRUE)
  pixelatorR:::assert_single_value(tidy, "bool")

  if (is.null(formula_str) && !is.null(vars)) {
    # Build the formula string
    formula_str <- paste(response, "~", paste(vars, collapse = "*"))
  } else if (is.null(formula_str) && is.null(vars)) {
    # If no vars are provided, use the response variable only
    formula_str <- paste(response, "~ 1")
  }

  # Run model
  mdl <- try(
    {
      aov(as.formula(formula_str), data = data)
    },
    silent = TRUE
  )
  if (inherits(mdl, "try-error")) {
    return(NULL)
  }

  # Tidy up
  if (tidy) {
    mdl <-
      mdl %>%
      tidy()

    # Add residuals if not present (no variance in data)
    if (!any(mdl$term == "Residuals")) {
      mdl <-
        mdl %>%
        bind_rows(tibble(
          term = "Residuals",
          df = NA,
          sumsq = 0,
          meansq = 0,
          statistic = NA,
          p.value = NA
        ))
    }

    mdl <-
      mdl %>%
      mutate(
        ss_total = sum(sumsq),
        ms_error = meansq[term == "Residuals"],
        eta_squared = sumsq / ss_total,
        omega_squared = (sumsq - df * ms_error) / (ss_total + ms_error),
        # Omega squared is not defined for residuals:
        omega_squared = ifelse(term == "Residuals", NA, omega_squared)
      ) %>%
      select(term, df, sumsq, meansq, statistic,
        p = p.value, eta_squared, omega_squared
      )
  }

  return(mdl)
}

#' Tidy variables by removing those with only one unique value
#'
#' This function selects variables from a data frame that have more than one unique value.
#'
#' @param df A data frame containing the variables to be checked.
#' @param vars A character vector of variable names to be checked.
#'
#' @return A character vector of variable names that have more than one unique value.
#'
#' @noRd
#'
tidy_vars <-
  function(df, vars) {
    select(df, all_of(vars)) %>%
      apply(2, n_distinct) %>%
      {
        which(. > 1)
      } %>%
      names()
  }

#' Run ANOVA for abundance data in a Seurat object
#'
#' This function runs an ANOVA for abundance data in a Seurat object, allowing for parallel processing.
#'
#' @param object A Seurat object containing the abundance data.
#' @param vars A character vector of length 1 or 2 specifying the variables to use in the ANOVA.
#' @param layer A string specifying the layer of the Seurat object to use for the analysis (default is "data").
#' @param mc_cores An integer specifying the number of cores to use for parallel processing (default is 1).
#' @param p_adj_method A string specifying the method for adjusting p-values (default is "bonferroni").
#'
#' @return A data frame with ANOVA results for each marker, including adjusted p-values.
#'
#' @export
#'
run_abundance_anova <-
  function(
    object,
    vars,
    layer = "data",
    mc_cores = 1,
    p_adj_method = "bonferroni"
  ) {
    pixelatorR:::assert_class(object, "Seurat")
    pixelatorR:::assert_vector(vars, "character", n = 1)
    pixelatorR:::assert_max_length(vars, 2)
    pixelatorR:::assert_single_value(layer, "string")
    pixelatorR:::assert_single_value(mc_cores, "numeric")
    pixelatorR:::assert_single_value(p_adj_method, "string")
    pixelatorR:::assert_is_one_of(p_adj_method, p.adjust.methods)

    cnt_data <-
      LayerData(object, layer) %>%
      as_tibble(rownames = "marker") %>%
      pivot_longer(
        cols = -marker,
        names_to = "comp_id",
        values_to = "x"
      )

    comp_meta_data <-
      FetchData(object, vars = vars) %>%
      as_tibble(rownames = "comp_id")

    vars <- tidy_vars(comp_meta_data, vars)
    if (length(vars) == 0) {
      cli::cli_abort("No valid variables provided for ANOVA. Variables with only one unique value were filtered out,
                     leaving no variables for comparison.")
    }

    aov_res <-
      cnt_data %>%
      left_join(comp_meta_data, by = "comp_id") %>%
      group_by(marker) %>%
      group_split() %>%
      parallel::mclapply(function(x) {
        .run_anova(x, vars, tidy = TRUE) %>%
          {
            if (is.null(.)) {
              .
            } else {
              mutate(., marker = x$marker[1]) %>%
                select(marker, everything())
            }
          }
      }, mc.cores = mc_cores) %>%
      bind_rows()

    aov_res <- aov_res %>%
      mutate(p_adj = p.adjust(p, method = p_adj_method)) %>%
      relocate(p_adj, .after = p)

    return(aov_res)
  }


#' Run ANOVA for proximity scores in a Seurat object
#'
#' This function runs an ANOVA for proximity scores in a Seurat object, allowing for parallel processing.
#'
#' @param proximity_scores A data frame containing proximity scores.
#' @param vars A character vector of length 1 or 2 specifying the variables to use in the ANOVA.
#' @param proximity_score A string specifying the column name for proximity scores (default is "log2_ratio").
#' @param mc_cores An integer specifying the number of cores to use for parallel processing (default is 1).
#' @param p_adj_method A string specifying the method for adjusting p-values (default is "bonferroni").
#'
#' @return A data frame with ANOVA results for each marker, including adjusted p-values.
#'
#' @export
#'
run_proximity_anova <-
  function(
    proximity_scores,
    vars,
    proximity_score = "log2_ratio",
    mc_cores = 1,
    p_adj_method = "bonferroni"
  ) {
    pixelatorR:::assert_class(proximity_scores, "tbl_df")
    pixelatorR:::assert_vector(vars, "character", n = 1)
    pixelatorR:::assert_max_length(vars, 2)
    pixelatorR:::assert_single_value(mc_cores, "numeric")
    pixelatorR:::assert_single_value(p_adj_method, "string")
    pixelatorR:::assert_is_one_of(p_adj_method, p.adjust.methods)

    for (var in vars) pixelatorR:::assert_col_in_data(var, proximity_scores)

    comp_meta_data <-
      proximity_scores %>%
      ungroup() %>%
      select(comp_id = sample_component, all_of(vars)) %>%
      distinct()

    vars <- tidy_vars(comp_meta_data, vars)
    if (length(vars) == 0) {
      cli::cli_abort("No valid variables provided for ANOVA. Variables with only one unique value were filtered out,
                     leaving no variables for comparison.")
    }

    proximity_scores_wide <-
      proximity_scores %>%
      rename(component = sample_component) %>%
      ProximityScoresToAssay(
        values_from = proximity_score,
        missing_obs = 0
      )

    # Ensure same order as metadata
    proximity_scores_wide <-
      proximity_scores_wide[, comp_meta_data$comp_id]

    aov_res <-
      seq_len(nrow(proximity_scores_wide)) %>%
      parallel::mclapply(function(i) {
        comp_meta_data %>%
          mutate(x = proximity_scores_wide[i, ]) %>%
          .run_anova(vars, tidy = TRUE) %>%
          {
            if (is.null(.)) {
              .
            } else {
              mutate(., contrast = rownames(proximity_scores_wide)[i]) %>%
                select(contrast, everything())
            }
          }
      }, mc.cores = mc_cores) %>%
      bind_rows()

    aov_res <- aov_res %>%
      separate(contrast, into = c("marker_1", "marker_2"), sep = "/") %>%
      mutate(p_adj = p.adjust(p, method = p_adj_method))

    return(aov_res)
  }


#' Find top proximity markers based on standard deviation of mean log2 ratio
#'
#' This function identifies the top markers based on the standard deviation of their mean log2 ratio
#' across samples, filtering by minimum percentage detected and range.
#'
#' @param proximity_score_summary A data frame containing proximity score mean summaries.
#' @param n_markers An integer specifying the number of top markers to select (default is 40).
#' @param min_pct_detected A numeric value specifying the minimum percentage of detection required (default is 0.25).
#' @param min_range A numeric value specifying the minimum range of mean log2 ratio required (default is 0.2).
#'
#' @return A character vector of the top marker names.
#'
#' @noRd
#'
find_top_proximity_markers <-
  function(
    proximity_score_summary,
    n_markers = 40,
    min_pct_detected = 0.25,
    min_range = 0.2
    ) {

    proximity_score_summary %>%
      group_by(sample_alias, marker_1) %>%
      # ensure all possible marker_2 per sample-marker_1 combo
      complete(marker_2 = unique(proximity_score_summary$marker_2)) %>%
      mutate(mean_log2_ratio = ifelse(pct_detected >= min_pct_detected, mean_log2_ratio, 0),
             mean_log2_ratio = ifelse(is.na(mean_log2_ratio), 0, mean_log2_ratio)) %>%
      group_by(marker_1) %>%
      summarize(
        sd = sd(mean_log2_ratio),
        max = max(mean_log2_ratio),
        min = min(mean_log2_ratio)
      ) %>%
      ungroup() %>%
      filter(max >= min_range | min <= -min_range) %>%
      top_n(n = n_markers, wt = sd) %>%
      slice_head(n = n_markers) %>%
      select(-sd) %>%
      pull(marker_1)
  }


#' Find top abundance markers based on mean expression
#'
#' This function identifies the top abundance markers based on mean abundance across samples and cell types.
#'
#' @param pg_data_processed A Seurat object containing processed data.
#' @param n_markers An integer specifying the number of top markers to select (default is 40).
#' @param summary_method A string specifying the summary method to use ("max" or "mean"; default is "max").
#' @param group_col An optional string specifying the column to group by (default is NULL).
#'
#' @return A data frame with the top markers for each cell type.
#'
#' @noRd
#'
find_top_abundance_markers <-
  function(
    pg_data_processed,
    n_markers = 40,
    summary_method = c("max", "mean"),
    group_col = NULL
  ) {
    summary_method <- match.arg(summary_method, c("max", "mean"))

    norm_data <-
      LayerData(pg_data_processed, "data")

    cell_annotation <-
      FetchData(pg_data_processed,
                vars = c("sample_alias", group_col)) %>%
      as_tibble(rownames = "cell_id")


    cell_annotation %>%
      group_by(across(all_of(c("sample_alias", group_col)))) %>%
      reframe(mean = enframe(rowMeans(norm_data[, cell_id, drop = FALSE]),
                             "marker", "mean")) %>%
      unnest(cols = c(mean)) %>%
      group_by(across(all_of(c("marker", group_col)))) %>%
      summarize(max = max(mean),
                mean = mean(mean),
                .groups = "drop") %>%
      group_by(across(all_of(c(group_col)))) %>%
      arrange(desc(!!sym(summary_method))) %>%
      slice_head(n = n_markers) %>%
      ungroup() %>%
      select(any_of(c("marker", group_col)))

  }
