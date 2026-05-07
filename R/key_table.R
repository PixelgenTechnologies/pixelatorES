#' Get top markers for each group in the object
#'
#' Calculate the top abundant markers for each group in the object and report the fraction of counts of the top markers
#' across all cells in the group.
#'
#' @param object A Seurat object.
#' @param group A character string specifying the group variable in the object.
#' @param n_markers A numeric vector specifying the number of top markers to report for each group.
#'
#' @return A tibble with the group ID, top markers, and their fractions.
#'
#' @export
get_top_markers <-
  function(object,
           group,
           n_markers = c(3, 5)) {
    pixelatorR:::assert_class(object, "Seurat")
    pixelatorR:::assert_single_value(group, type = "string")
    pixelatorR:::assert_vector(n_markers, "numeric", n = 1)

    sample_cells <-
      FetchData(object, group) %>%
      as_tibble(rownames = "cell_id") %>%
      rename(group_id = !!sym(group)) %>%
      group_by(group_id)

    cts <- LayerData(object, "counts")


    sample_cells %>%
      group_split() %>%
      set_names(group_keys(sample_cells)$group_id) %>%
      lapply(function(cells) {
        prop_table <-
          cts[, cells$cell_id] %>%
          as.matrix() %>%
          rowSums() %>%
          prop.table() %>%
          enframe() %>%
          arrange(desc(value))

        n_markers %>%
          set_names(paste0("top", ., "_fraction")) %>%
          lapply(function(n) {
            sum(prop_table$value[1:n])
          }) %>%
          as_tibble() %>%
          mutate(top_markers = paste(prop_table$name[1:max(n_markers)], collapse = ", "))
      }) %>%
      bind_rows(.id = group)
  }

#' Get sequencing saturation metrics
#'
#' This function calculates sequencing saturation metrics for each sample based on the provided quality control metrics.
#'
#' @param object A Seurat object.
#' @param sample_qc_metrics A list containing quality control metrics for each sample.
#'
#' @return A tibble containing sequencing saturation metrics for each sample.
#'
#' @export
#'
get_seq_saturation <-
  function(object, sample_qc_metrics) {
    pixelatorR:::assert_class(sample_qc_metrics, "list")

    if (!is.null(sample_qc_metrics$pool_qc_files)) {
      group_col <- "pool"

      fetched_data <-
        sample_qc_metrics$pool_qc_files %>%
        lapply(function(sample_qc_data) {
          tibble(
            total_reads = sample_qc_data$amplicon$input_reads,
            deduped_valid_reads = sample_qc_data$collapse$output_molecules,
            valid_reads = sample_qc_data$collapse$input_reads
          )
        })
    } else {
      group_col <- "sample_alias"

      fetched_data <-
        sample_qc_metrics$qc_files %>%
        lapply(function(sample_qc_data) {
          tibble(
            Q30 = sample_qc_data$amplicon$q30_statistics$total,
            total_reads = sample_qc_data$amplicon$input_reads,
            deduped_valid_reads = sample_qc_data$collapse$output_molecules,
            valid_reads = sample_qc_data$collapse$input_reads
          )
        })
    }

    summarised_object_counts <-
      object[[]] %>%
      as_tibble() %>%
      group_by(!!sym(group_col)) %>%
      summarise(
        graph_edges = sum(n_edges),
        graph_proteins = sum(n_umi),
        graph_reads = sum(reads_in_component)
      )

    fetched_data %>%
      bind_rows(.id = group_col) %>%
      inner_join(
        summarised_object_counts,
        by = group_col
      ) %>%
      mutate(
        valid_reads_saturation = sequencing_saturation(
          deduped_valid_reads,
          valid_reads
        ),
        graph_edge_saturation = sequencing_saturation(
          graph_edges,
          graph_reads
        ),
        graph_node_saturation = sequencing_saturation(
          graph_proteins,
          graph_reads
        ),
        fraction_valid_reads = 100 * valid_reads / total_reads,
        fraction_graph_reads = 100 * graph_reads / total_reads
      )
  }


#' Get crossing edges data
#'
#' This function retrieves the crossing edges data from the sample quality control metrics.
#'
#' @param sample_qc_metrics A list containing quality control metrics for each sample.
#'
#' @return A tibble containing the crossing edges data for each sample.
#'
#' @export
#'
get_crossing_edges <-
  function(sample_qc_metrics) {
    pixelatorR:::assert_class(sample_qc_metrics, "list")

    if (!is.null(sample_qc_metrics$pool_qc_files)) {
      group_col <- "pool"

      fetched_data <-
        sample_qc_metrics$pool_qc_files %>%
        lapply(function(sample_qc_data) {
          tibble(
            `Initial stage` =
              sample_qc_data$graph$crossing_edges_removed_initial_stage,
            `Refinement stage` = sample_qc_data$graph$crossing_edges_removed -
              sample_qc_data$graph$crossing_edges_removed_initial_stage,
            removed_total = sample_qc_data$graph$crossing_edges_removed,
            total_edges_in = sample_qc_data$graph$edge_count_pre_recovery
          )
        })
    } else {
      group_col <- "sample_alias"

      fetched_data <-
        sample_qc_metrics$qc_files %>%
        lapply(function(sample_qc_data) {
          tibble(
            `Initial stage` =
              sample_qc_data$graph$crossing_edges_removed_initial_stage,
            `Refinement stage` = sample_qc_data$graph$crossing_edges_removed -
              sample_qc_data$graph$crossing_edges_removed_initial_stage,
            removed_total = sample_qc_data$graph$crossing_edges_removed,
            total_edges_in = sample_qc_data$graph$edge_count_pre_recovery
          )
        })
    }

    fetched_data %>%
      bind_rows(.id = group_col) %>%
      pivot_longer(
        cols = c("Initial stage", "Refinement stage"),
        names_to = "type", values_to = "edges"
      ) %>%
      mutate(
        percent = edges / total_edges_in * 100
      )
  }

#' Get denoising data
#'
#' This function retrieves the denoising data from the sample quality control metrics.
#'
#' @param sample_qc_metrics A list containing quality control metrics for each sample.
#'
#' @return A tibble containing the denoising data for each sample.
#'
#' @export
#'
get_denoising_data <-
  function(sample_qc_metrics) {
    pixelatorR:::assert_class(sample_qc_metrics, "list")

    if (!is.null(sample_qc_metrics$qc_files)) {
      sample_qc_metrics <- sample_qc_metrics$qc_files
    }

    lapply(names(sample_qc_metrics), function(nm) {
      tibble(
        sample_alias = nm,
        ratio =
          sample_qc_metrics[[nm]]$denoise$ratio_of_umis_removed * 100,
        number_of_umis_removed =
          sample_qc_metrics[[nm]]$denoise$number_of_umis_removed
      )
    }) %>%
      bind_rows()
  }

#' Get coreness data
#'
#' This function retrieves coreness data from a Seurat object and summarizes it by sample and component.
#'
#' @param object A Seurat object containing the coreness data.'
#'
#' @return A list containing the coreness data, component summary, and sample summary.
#'
#' @export
#'
get_coreness_data <-
  function(object) {
    pixelatorR:::assert_class(object, "Seurat")

    # Get k coreness data
    k_coreness_data <-
      object[[]] %>%
      as_tibble(rownames = "sample_component") %>%
      select(1, sample_alias, molecules = n_umi, matches("^k_core")) %>%
      pivot_longer(
        cols = matches("k_core"),
        names_to = "k_core", values_to = "nodes"
      ) %>%
      mutate(
        nodes = ifelse(!is.finite(nodes), 0, nodes),
        percent_nodes = 100 * nodes / molecules,
        coreness = as.numeric(str_extract(k_core, "[0-9]+"))
      )

    k_coreness_summary <-
      k_coreness_data %>%
      group_by(sample_component, sample_alias) %>%
      summarise(
        mean_coreness = weighted.mean(coreness, nodes),
        percent_dangling_nodes = percent_nodes[coreness == 1],
        percent_well_connected_nodes = sum(percent_nodes[coreness %in% 3:5])
      )

    k_coreness_summary_sample <-
      k_coreness_summary %>%
      group_by(sample_alias) %>%
      summarise(
        median_mean_coreness = median(mean_coreness),
        median_percent_dangling_nodes = median(percent_dangling_nodes),
        median_percent_well_connected_nodes = median(percent_well_connected_nodes)
      )

    return_list <-
      list(
        data = k_coreness_data,
        component_summary = k_coreness_summary,
        sample_summary = k_coreness_summary_sample
      )

    return(return_list)
  }

#' Get degree distribution data
#'
#' This function retrieves degree distribution data from the sample quality control metrics.
#'
#' @param sample_qc_metrics A list containing quality control metrics for each sample.
#'
#' @return A tibble containing the degree distribution data for each sample.
#'
#' @export
#'
get_degree_distribution <-
  function(sample_qc_metrics) {
    pixelatorR:::assert_class(sample_qc_metrics, "list")

    if (!is.null(sample_qc_metrics$pool_qc_files)) {
      group_col <- "pool"

      qc_data <-
        sample_qc_metrics$pool_qc_files
    } else {
      group_col <- "sample_alias"

      qc_data <-
        sample_qc_metrics$qc_files
    }

    qc_data %>%
      lapply(function(sample_qc_data) {
        list(
          umi1 = sample_qc_data$collapse$umi1_degree_distribution,
          umi2 = sample_qc_data$collapse$umi2_degree_distribution
        ) %>%
          map(
            . %>%
              unlist() %>%
              enframe("degree", "n")
          ) %>%
          bind_rows(.id = "umi_type")
      }) %>%
      bind_rows(.id = group_col) %>%
      group_by(!!sym(group_col)) %>%
      mutate(
        degree = as.integer(degree),
        percent_nodes = 100 * n / sum(n)
      ) %>%
      ungroup()
  }

#' Get read statistics for each sample
#'
#' This function retrieves read statistics for each sample from a Seurat object.
#'
#' @param object A Seurat object containing the data.
#'
#' @return A tibble summarizing read statistics for each sample.
#'
#' @export
#'
get_read_stats <-
  function(object) {
    pixelatorR:::assert_class(object, "Seurat")

    FetchData(object, c(
      "sample_alias", "reads_in_component", "n_umi",
      "isotype_fraction", "intracellular_fraction"
    )) %>%
      group_by(sample_alias) %>%
      summarise(
        # Number of cells
        n_cells = n(),
        n_cells_over10k = sum(n_umi > 10000),

        # Reads per cell
        median_reads_per_cell = median(reads_in_component),

        # Abs per cell
        median_abs_per_cell = median(n_umi),

        # Isotype and control fraction
        median_isotype_count_pct = 100 * median(isotype_fraction),
        median_intracellular_count_pct = 100 * median(intracellular_fraction),
      )
  }


#' Get sample hashing statistics for each sample
#'
#' @param object A Seurat object containing the data.
#' @param sample_qc_metrics A list containing quality control metrics for each sample, which may include hashing
#' information.
#'
#' @return A list of summary tables, or `NULL` when no hash count columns exist in metadata.
#'
#' @export
#'
get_hash_stats <- function(object, sample_qc_metrics) {
  pixelatorR:::assert_class(object, "Seurat")

  hash_counts <- object[[]] %>%
    select(matches("hash_counts"))

  if (ncol(hash_counts) == 0) {
    return(NULL)
  }

  sample_hash_counts <- tibble(
    sample_alias = object$sample_alias,
    hash_counts = rowSums(hash_counts),
    total_counts = object$n_umi
  ) %>%
    group_by(sample_alias) %>%
    summarise(
      total_hash_counts = sum(hash_counts),
      total_counts = sum(total_counts),
      .groups = "drop"
    ) %>%
    mutate(hash_pct = total_hash_counts / total_counts)

  id_part <- str_extract(colnames(hash_counts)[seq_len(ncol(hash_counts))], "(?<=_)[[:alnum:]]+(?=\\.)")
  version_part <- str_extract(colnames(hash_counts)[seq_len(ncol(hash_counts))], "(?<=\\.)\\d+$") %>% as.integer()

  component_hash_counts <- lapply(unique(id_part), function(id) {
    hash_counts_subset <- hash_counts[, id_part == id, drop = FALSE]
    hash_counts_subset <- hash_counts_subset[, order(version_part[id_part == id]), drop = FALSE]
    colnames(hash_counts_subset) <- paste0(id, "-", seq_len(ncol(hash_counts_subset)))
    hash_counts_subset <- cbind(hash_counts_subset, sample_alias = object$sample_alias) %>%
      rownames_to_column("component") %>%
      pivot_longer(cols = -all_of(c("component", "sample_alias")), names_to = "hash", values_to = "count") %>%
      separate(hash, into = c("id", "version"), sep = "-")
    return(hash_counts_subset)
  }) %>%
    bind_rows()

  component_stats <- component_hash_counts %>%
    group_by(component, sample_alias, id) %>%
    mutate(purity = count / sum(count)) %>%
    mutate(purity = ifelse(is.na(purity), 0, purity))

  n_ids <- id_part %>%
    unique() %>%
    length()
  n_hashes <- version_part %>%
    unique() %>%
    length()
  id_hash_order <- paste0(rep(unique(id_part), n_hashes), "-", rep(1:n_hashes, each = n_ids))

  .create_heatmap_df <- function(component_stats, var = "purity", sample_levels) {
    component_stats %>%
      unite("hash", id, version, sep = "-") %>%
      select(component, sample_alias, hash, !!sym(var)) %>%
      pivot_wider(names_from = hash, values_from = !!sym(var)) %>%
      ungroup() %>%
      mutate(sample_alias = factor(sample_alias, levels = sample_levels)) %>%
      rename(sample_component = component)
  }


  component_stats_heatmap_purity <- .create_heatmap_df(
    component_stats,
    "purity",
    levels(object$sample_alias)
  ) %>%
    left_join(
      object[[]] %>%
        select(n_umi) %>%
        rownames_to_column("sample_component"),
      by = "sample_component"
    ) %>%
    rename(molecules = n_umi)
  attr(component_stats_heatmap_purity, "id_hash_order") <- id_hash_order

  component_stats_heatmap_fraction <- .create_heatmap_df(
    component_stats,
    "count",
    levels(object$sample_alias)
  )

  numeric_cols <- names(component_stats_heatmap_fraction)[
    sapply(component_stats_heatmap_fraction, is.numeric)
  ]

  # Heatmap data with fractions of counts per component
  component_stats_heatmap_fraction <- component_stats_heatmap_fraction %>%
    rowwise() %>%
    mutate(
      r_sum = sum(dplyr::c_across(tidyselect::all_of(numeric_cols)), na.rm = TRUE)
    ) %>%
    mutate(dplyr::across(tidyselect::all_of(numeric_cols), ~ .x / r_sum)) %>%
    select(-r_sum)
  attr(component_stats_heatmap_fraction, "id_hash_order") <- id_hash_order

  sample_stats <- component_stats %>%
    group_by(component, id) %>%
    dplyr::slice_max(order_by = purity, n = 1, with_ties = FALSE) %>%
    group_by(sample_alias, id) %>%
    summarise(
      mean_purity = mean(purity),
      .groups = "drop"
    ) %>%
    mutate(sample_alias = factor(sample_alias, levels = unique(object$sample_alias))) %>%
    pivot_wider(names_from = id, values_from = mean_purity, names_prefix = "mean_purity_") %>%
    left_join(sample_hash_counts %>% select(sample_alias, hash_pct), by = "sample_alias")

  # Sample confidence per component and sample
  component_sample_confidence <-
    sample_qc_metrics$pool_qc_files %>%
    lapply(function(pool) {
      pool$sample_calling$sample_confidences_per_sample %>%
        lapply(tibble) %>%
        bind_rows(.id = "sample_alias") %>%
        rename(sample_confidence = 2)
    }) %>%
    bind_rows(.id = "pool") %>%
    mutate(
      sample_alias =
        ifelse(str_detect(sample_alias, "undetermined$"),
          str_remove(sample_alias, paste0(pool, "_")),
          sample_alias
        )
    )

  list(
    component_stats = component_stats,
    component_stats_heatmap_purity = component_stats_heatmap_purity,
    component_stats_heatmap_fraction = component_stats_heatmap_fraction,
    sample_stats = sample_stats,
    component_sample_confidence = component_sample_confidence
  )
}


#' Get quality control metrics for samples
#'
#' This function retrieves quality control metrics for samples from the provided sample QC metrics and sample sheet.
#'
#' @param object A Seurat object containing the data.
#' @param sample_qc_metrics A list containing quality control metrics for each sample.
#' @param sample_sheet A data frame containing sample metadata, including `sample_alias`.
#'
#' @return A list containing various quality control metrics for each sample.
#'
#' @export
#'
get_qc_metrics <-
  function(object, sample_qc_metrics, sample_sheet) {
    .format <-
      function(tb) {
        if (is.null(tb)) {
          return(NULL)
        }
        if ("sample_alias" %in% names(tb)) {
          sample_levels <-
            sample_sheet$sample_alias[sample_sheet$sample_alias %in% unique(tb$sample_alias)]
          if ("undetermined" %in% tb$sample_alias) sample_levels <- c(sample_levels, "undetermined")
          tb <- order_sample_alias_factors(tb, levels = sample_levels)
        }
        if ("pool" %in% names(tb)) {
          pool_levels <-
            unique(sample_sheet$pool[!is.na(sample_sheet$pool)])
          pool_levels <- pool_levels[pool_levels %in% unique(tb$pool)]
          tb <- order_sample_alias_factors(tb,
            levels = pool_levels,
            column_name = "pool"
          )
        }
        if (inherits(tb, "list")) {
          tb <- lapply(tb, .format)
        }

        return(tb)
      }

    list(
      read_stats = get_read_stats(object),
      sample_hash_stats = get_hash_stats(object, sample_qc_metrics),
      seq_saturation = get_seq_saturation(object, sample_qc_metrics),
      crossing_edges = get_crossing_edges(sample_qc_metrics),
      degree_distribution = get_degree_distribution(sample_qc_metrics),
      denoising = get_denoising_data(sample_qc_metrics),
      coreness = get_coreness_data(object),
      top_markers = get_top_markers(object, "sample_alias",
        n_markers = c(3, 5)
      )
    ) %>%
      lapply(.format)
  }

#' Process metrics table content
#'
#' @param table_content_list List of metric tables to join.
#' @param detailed Whether to restrict to a short set of columns.
#'
#' @noRd
#'
.process_metrics_table_content <-
  function(table_content_list, detailed = TRUE) {
    join_by <- ifelse("sample_alias" %in% names(table_content_list[[1]]),
      "sample_alias", "pool"
    )

    table_content <-
      table_content_list %>%
      keep(~ !is.null(.x)) %>%
      reduce(left_join, by = join_by) %>%
      select(1, any_of(key_metric_definitions$var))

    if (!detailed) {
      table_content <-
        table_content %>%
        select(1, any_of(c(
          "n_cells",
          "median_isotype_count_pct",
          "median_abs_per_cell",
          "median_reads_per_cell",
          "total_reads",
          "graph_node_saturation",
          "graph_edge_saturation",
          "median_mean_coreness"
        )))
    }

    for (i in seq(from = 2, to = ncol(table_content))) {
      definitions_i <- which(key_metric_definitions$var == colnames(table_content)[i])

      if (is.numeric(table_content[[i]])) {
        table_content[, i] <- table_content[, i] / key_metric_definitions$scale[definitions_i]
      }
      colnames(table_content)[i] <- key_metric_definitions$display_name[definitions_i]
    }

    table_content %>%
      mutate(across(where(is.numeric), . %>%
        round(2))) %>%
      mutate_all(as.character)
  }


#' Create a key metric table for samples
#'
#' This function creates a key metric table for samples by combining various
#' quality control metrics and formatting them for display.
#'
#' @param qc_metrics_tables A list containing quality control metrics tables
#' for each sample.
#' @param detailed A logical value indicating whether to include all metrics.
#' @param return_data A logical value indicating whether to return the data instead of the formatted table.
#'
#' @return A named list with elements `sample` and `pool` (the latter may be `NULL`). With
#'   `return_data = TRUE`, returns pre-styled tibbles in the same structure.
#'
#' @export
#'
key_metric_table <-
  function(qc_metrics_tables, detailed = TRUE, return_data = FALSE) {
    table_content_list <-
      list(
        qc_metrics_tables$read_stats,
        qc_metrics_tables$seq_saturation,
        qc_metrics_tables$denoising,
        qc_metrics_tables$coreness$sample_summary,
        qc_metrics_tables$top_markers,
        qc_metrics_tables$crossing_edges %>%
          select(
            any_of(c("sample_alias", "pool")),
            type,
            percent
          ) %>%
          pivot_wider(names_from = "type", values_from = "percent")
      )

    if ("sample_hash_stats" %in% names(qc_metrics_tables) && !is.null(qc_metrics_tables$sample_hash_stats)) {
      table_content_list <- c(table_content_list, list(qc_metrics_tables$sample_hash_stats$sample_stats))
    }

    pool_content <-
      sapply(table_content_list, function(x) !is.null(x) && "pool" %in% names(x))

    sample_content <-
      !pool_content

    table_content <-
      list(
        sample = table_content_list[sample_content],
        pool = table_content_list[pool_content]
      )
    table_content <- table_content[which(sapply(table_content, function(x) length(x) != 0))]
    table_content <- lapply(table_content, .process_metrics_table_content, detailed = detailed)

    table_content$sample <-
      table_content$sample %>%
      rename(
        `Sample ID` = sample_alias
      )

    if (!is.null(table_content$pool)) {
      table_content$pool <-
        table_content$pool %>%
        rename(
          `Pool ID` = pool
        )
    }

    if (return_data) {
      return(table_content)
    }

    table_content <-
      table_content %>%
      purrr::map(. %>%
        pivot_longer(
          cols = -1,
          names_to = "Metric",
          values_to = "Value"
        ) %>%
        pivot_wider(
          names_from = 1,
          values_from = "Value"
        ) %>%
        left_join(
          key_metric_definitions %>%
            select(display_name, description),
          by = c("Metric" = "display_name")
        ) %>%
        mutate(
          Metric = format_with_info_bootstrap(Metric, description)
        ) %>%
        select(-description))

    if (!is.null(table_content$pool)) {
      pool_table <-
        table_content$pool %>%
        style_table(
          escape = FALSE,
          tooltips = TRUE
        )
    } else {
      pool_table <- NULL
    }

    sample_table <-
      table_content$sample %>%
      style_table(
        escape = FALSE,
        tooltips = TRUE
      )

    list(
      pool = pool_table,
      sample = sample_table
    )
  }


#' Format content with an info icon and tooltip
#'
#' This function formats content with an info icon that displays a tooltip
#' with additional information when hovered over.
#'
#' @param content A string containing the content to be formatted.
#' @param description A string containing the description to be displayed in the tooltip.
#'
#' @return A formatted string with an info icon and tooltip.
#'
#' @export
#'
format_with_info_bootstrap <-
  function(content, description) {
    # Use Bootstrap info icon with tooltip
    sprintf(
      '%s <i class="bi bi-info-circle-fill text-primary"
        data-bs-toggle="tooltip"
        data-bs-title="%s"
        style="cursor:help;"></i>',
      content,
      description
    )
  }

#' Key metric definitions
#'
#' This tibble contains definitions for key metrics.
#'
#' @noRd
#'
key_metric_definitions <-
  tibble::tribble(
    ~var, ~display_name, ~scale,
    ~description,
    "n_cells", "Number of cells", 1,
    "Total number of cells in the sample.",
    "n_cells_over10k", "Number of cells >10k nodes", 1,
    "Number of cells with more than 10,000 nodes (proteins).",
    "median_isotype_count_pct", "Median isotype % counts", 1,
    "Median percentage of counts attributed to isotype controls across all cells in the sample.",
    "median_intracellular_count_pct", "Median intracellular % counts", 1,
    "Median percentage of counts attributed to intracellular proteins across all cells in the sample.",
    "median_abs_per_cell", "Median proteins per cell [k]", 1e3,
    "Median number of proteins (nodes) per cell, scaled to thousands.",
    "median_reads_per_cell", "Median reads per cell [k]", 1e3,
    "Median number of reads per cell, scaled to thousands.",
    "Q30", "Q30 [%]", 1e-2,
    "Percentage of bases with a Q30 Phred score of at least 30. Q30 indicates a 99.9% base call accuracy.",
    "total_reads", "Total reads [M]", 1e6,
    "Total number of reads in the sample, scaled to millions.",
    "valid_reads", "Valid reads [M]", 1e6,
    "Total number of valid reads in the sample, excluding reads that do not map to any panel protein,
    scaled to millions.",
    "graph_proteins", "Graph Nodes [M]", 1e6,
    "Total number of graph nodes (proteins) in the sample, excluding nodes that don't belong to any cell,
    scaled to millions.",
    "graph_edges", "Graph Edges [M]", 1e6,
    "Total number of graph edges in the sample, excluding edges that don't belong to any cell, scaled to millions.",
    "graph_node_saturation", "Graph node saturation [%]", 1,
    "Percentage of graph nodes that are saturated, indicating the proportion of nodes that have been fully sampled.",
    "graph_edge_saturation", "Graph edge saturation [%]", 1,
    "Percentage of graph edges that are saturated, indicating the proportion of edges that have been fully sampled.",
    "valid_reads_saturation", "Valid reads saturation [%]", 1,
    "Percentage of valid reads that are saturated, indicating the proportion of reads that have been fully sampled.",
    "fraction_valid_reads", "Valid reads fraction [%]", 1,
    "Percentage of total reads that are valid reads.",
    "fraction_graph_reads", "Graph reads fraction [%]", 1,
    "Percentage of total reads that are graph reads,
    indicating the proportion of reads that contribute to cell graphs.",
    "ratio", "% Denoised UMIs", 1,
    "Percentage of unique molecular identifiers (UMIs) that have been removed in denoising.",
    "number_of_umis_removed", "Total denoised UMIs [M]", 1e6,
    "Total number of unique molecular identifiers (UMIs) that have been removed during the denoising process,
    scaled to millions.",
    "median_mean_coreness", "Median mean coreness", 1,
    "Median mean coreness across all components in the sample,
    indicating the average connectivity of nodes in cell graphs.",
    "median_percent_dangling_nodes", "Median % dangling nodes", 1,
    "Median percentage of dangling nodes in the sample,
    indicating the proportion of nodes that have a coreness of 1.",
    "median_percent_well_connected_nodes", "Median % well connected nodes", 1,
    "Median percentage of well-connected nodes in the sample,
    indicating the proportion of nodes that have a coreness of 3-5.",
    "top3_fraction", "Top 3 % counts", 1e-2,
    "Percentage of counts attributed to the top 3 abundant markers in the sample.",
    "top5_fraction", "Top 5 % counts", 1e-2,
    "Percentage of counts attributed to the top 5 abundant markers in the sample.",
    "top_markers", "Top 5 markers", 1,
    "The top 5 abundant markers in the sample.",
    "Initial stage", "% Crossing edges (Initial)", 1,
    "Percentage of edges that were crossing edges in the initial stage of graph refinement.",
    "Refinement stage", "% Crossing edges (Refinement)", 1,
    "Percentage of edges that were crossing edges in the refinement stage of graph refinement.",
    "mean_purity_B2M", "% B2M hash purity", 1e-2,
    "Mean purity of the B2M hash across components in the sample.",
    "mean_purity_CD298", "% CD298 hash purity", 1e-2,
    "Mean purity of the CD298 hash across components in the sample.",
    "mean_purity_CD98", "% CD98 hash purity", 1e-2,
    "Mean purity of the CD98 hash across components in the sample.",
    "hash_pct", "% hash counts", 1e-2,
    "Percentage of total UMI counts that are attributed to hashing antibodies."
  )
