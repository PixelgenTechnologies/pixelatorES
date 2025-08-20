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

    sample_qc_metrics %>%
      lapply(function(sample_qc_data) {
        tibble(
          total_reads = sample_qc_data$amplicon$input_reads,
          deduped_valid_reads = sample_qc_data$collapse$output_molecules,
          valid_reads = sample_qc_data$collapse$input_reads
        )
      }) %>%
      bind_rows(.id = "sample_alias") %>%
      left_join(
        object[[]] %>%
          as_tibble() %>%
          group_by(sample_alias) %>%
          summarise(
            graph_edges = sum(n_edges),
            graph_proteins = sum(n_umi),
            graph_reads = sum(reads_in_component)
          ),
        by = c("sample_alias")
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

    sample_qc_metrics %>%
      lapply(function(sample_qc_data) {
        tibble(
          `Initial stage` =
            sample_qc_data$graph$crossing_edges_removed_initial_stage,
          `Refinement stage` = sample_qc_data$graph$crossing_edges_removed -
            sample_qc_data$graph$crossing_edges_removed_initial_stage,
          removed_total = sample_qc_data$graph$crossing_edges_removed,
          total_edges_in = sample_qc_data$graph$edge_count_pre_recovery
        )
      }) %>%
      bind_rows(.id = "sample_alias") %>%
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

    # Get denoising data
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
    list(
      read_stats = get_read_stats(object),
      seq_saturation = get_seq_saturation(object, sample_qc_metrics),
      crossing_edges = get_crossing_edges(sample_qc_metrics),
      denoising = get_denoising_data(sample_qc_metrics),
      coreness = get_coreness_data(object),
      top_markers = get_top_markers(object, "sample_alias",
        n_markers = c(3, 5)
      )
    ) %>%
      lapply(order_sample_alias_factors, levels = sample_sheet$sample_alias)
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
#' @return A formatted table of key metrics for each sample.
#'
#' @export
#'
key_metric_table <-
  function(qc_metrics_tables, detailed = TRUE, return_data = FALSE) {
    table_content <-
      list(
        # Basic cell stats
        qc_metrics_tables$read_stats,

        # Sequencing saturation
        qc_metrics_tables$seq_saturation,

        # Denoise stats
        qc_metrics_tables$denoising,

        # Coreness stats
        qc_metrics_tables$coreness$sample_summary,

        # Top markers
        qc_metrics_tables$top_markers,

        # Crossing edges
        qc_metrics_tables$crossing_edges %>%
          select(
            sample_alias,
            type,
            percent
          ) %>%
          pivot_wider(names_from = "type", values_from = "percent")
      ) %>%
      reduce(left_join, by = "sample_alias") %>%
      select(1, any_of(key_metric_definitions$var))

    if (!detailed) {
      table_content <-
        table_content %>%
        select(1, all_of(c(
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


    table_content <-
      table_content %>%
      # Format
      mutate(across(where(is.numeric), . %>%
        round(2))) %>%
      mutate_all(as.character) %>%
      rename(
        `Sample ID` = sample_alias
      )


    if (return_data) {
      return(table_content)
    }

    table_content %>%
      # Pivot
      pivot_longer(
        cols = -`Sample ID`,
        names_to = "Metric",
        values_to = "Value"
      ) %>%
      pivot_wider(
        names_from = "Sample ID",
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
      select(-description) %>%
      # Show table
      style_table(
        escape = FALSE,
        tooltips = TRUE
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
    "Percentage of edges that were crossing edges in the refinement stage of graph refinement."
  )
