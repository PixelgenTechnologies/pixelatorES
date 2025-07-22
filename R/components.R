#' Create the component for bleedover noise
#'
#' This function creates a ggplot object that visualizes the fraction of UMIs
#' removed by denoising.
#'
#' @param qc_metrics_tables A list containing a data frame named `denoising`
#'
#' @return A ggplot object visualizing the fraction of UMIs removed by denoising.
#'
#' @export
#'
component_bleedover_noise <- function(
  qc_metrics_tables) {
  pixelatorR:::assert_class(qc_metrics_tables, "list")
  p <-
    qc_metrics_tables$denoising %>%
    ggplot(aes(sample_alias, ratio)) +
    geom_col(fill = "#DAD6D7") +
    geom_text(aes(label = paste0(round(ratio, 3), " %")),
      vjust = -.1,
      size = 3
    ) +
    scale_y_continuous(
      expand = expansion(c(0, 0.15)),
      limits = c(0, 100),
      breaks = c(0, 25, 50, 75, 100)
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 60, hjust = 1),
      panel.grid = element_blank()
    ) +
    labs(
      x = NULL, y = "Removed UMIs [%]",
      title = "Fraction of UMIs removed by bleedover denoising",
      caption = expression(~ italic(frac("# Discarded UMIs", "# Total UMIs")))
    )
  return(p)
}

#' Create the component for control markers
#'
#' This function creates a list with two plots and a summary table.
#'
#' @param pg_data A Seurat object containing the data to be plotted.
#'
#' @return A list containing two violin plots and a summary table.
#'
#' @export
#'
component_control_markers <- function(
  pg_data) {
  pixelatorR:::assert_class(pg_data, "Seurat")
  plot_data <-
    pg_data[[]] %>%
    as_tibble() %>%
    mutate(
      isotype_counts = round(isotype_fraction * n_umi, 0) + 1
    )

  p1 <-
    plot_data %>%
    plot_violin(
      x = "sample_alias",
      y = "isotype_fraction",
      y_label = "Percent control markers",
      round = 4,
      expand = c(0, 0.2),
      use_pct = TRUE
    )

  p2 <-
    plot_data %>%
    plot_violin(
      x = "sample_alias",
      y = "isotype_counts",
      y_label = "Control markers counts",
      round = 4,
      expand = c(0, 0.2),
      use_log10 = TRUE
    )

  tabl <-
    plot_data %>%
    group_by(sample_alias) %>%
    summarize(
      median_isotype_percent = median(isotype_fraction) * 100,
      median_isotype_counts = median(isotype_counts)
    ) %>%
    select(
      `Sample ID` = sample_alias,
      `Median percent control markers` = median_isotype_percent,
      `Median control markers counts` = median_isotype_counts
    ) %>%
    mutate(`Median percent control markers` = round(`Median percent control markers`, 3)) %>%
    style_table(caption = "Median percent control markers", interactive = FALSE)


  return(list(p1 = p1, p2 = p2, tabl = tabl))
}

#' Create the component for molecule rank plots
#'
#' This function generates a list of ggplot objects, each representing the
#' molecule rank plot for a specific sample in the provided Seurat object.
#'
#' @param pg_data A Seurat object containing the data to be plotted.
#' @param params A list of parameters, including `molecule_rank_cutoff`, which
#' is used to draw a horizontal line in the plots.
#'
#' @return A list of ggplot objects, each corresponding to a sample's molecule
#' rank plot.
#'
#' @export
#'
component_molecule_rank_plot <- function(
  pg_data,
  params) {
  plots <-
    FetchData(pg_data, "sample_alias") %>%
    pull(sample_alias) %>%
    levels() %>%
    set_names() %>%
    lapply(function(x) {
      pg_data %>%
        FetchData(c("sample_alias", "n_umi", "rank")) %>%
        arrange(sample_alias == x) %>%
        ggplot(aes(rank, n_umi, color = sample_alias == x)) +
        geom_point(size = 0.5, show.legend = FALSE) +
        geom_hline(
          yintercept = params$molecule_rank_cutoff, linetype = "dashed",
          color = "#E05573"
        ) +
        scale_x_log10() +
        scale_y_log10() +
        scale_color_manual(values = c("TRUE" = "black", "FALSE" = "gray80")) +
        theme_bw() +
        labs(
          x = "Component rank (by number of molecules)",
          y = "Number of molecules",
          title = "Molecule rank plot"
        )
    })

  return(plots)
}

#' Create the component for molecule plot
#'
#' This function generates a violin plot showing the distribution of the number
#' of molecules per sample, with a horizontal line indicating the molecule count cutoff.
#'
#' @param pg_data A Seurat object containing the data to be plotted.
#' @param params A list of parameters, including `molecule_rank_cutoff`, which
#' is used to draw a horizontal line in the plot.
#' @param sample_palette A named vector of colors for the samples.
#'
#' @return A ggplot object representing the violin plot of molecule counts.
#'
#' @export
#'
component_molecule_plot <- function(
  pg_data,
  params,
  sample_palette) {
  p <- pg_data[[]] %>%
    plot_violin(
      x = "sample_alias",
      y = "n_umi",
      fill = "condition",
      y_label = "Number of molecules",
      round = 0,
      use_log10 = TRUE,
      palette = sample_palette,
      alpha = 1
    ) +
    geom_hline(yintercept = params$molecule_rank_cutoff, linetype = "dashed", color = "#E05573")

  return(p)
}


#' Create the component sequencing reads and molecules component
#'
#' This function creates a component that visualizes the sequencing reads
#' and molecules processed at each stage of the sequencing pipeline.
#'
#' @param sample_qc_metrics A list of sample QC metrics.
#' @param sample_levels Optional vector of sample levels to order the samples in the plots.
#'
#' @return A list containing plots and a table.
#'
#' @export
#'
component_sequencing_reads_and_molecules <-
  function(sample_qc_metrics, sample_levels = NULL) {
    plot_data <-
      sample_qc_metrics %>%
      lapply(function(sample_qc_data) {
        list(
          "amplicon" = list(
            input = sample_qc_data$amplicon$input_reads,
            output = sample_qc_data$amplicon$output_reads
          ),
          "demux" = list(
            input = sample_qc_data$demux$input_reads,
            output = sample_qc_data$demux$output_reads
          ),
          "collapse" = list(
            input = sample_qc_data$collapse$input_reads,
            output = sample_qc_data$collapse$output_molecules
          ),
          "graph" = list(
            input = sample_qc_data$graph$molecules_input,
            output = sample_qc_data$graph$molecules_output
          )
        ) %>%
          map(
            . %>%
              map(
                . %>%
                  enframe()
              ) %>%
              bind_rows(.id = "type")
          ) %>%
          bind_rows(.id = "stage") %>%
          pivot_wider(names_from = "type", values_from = "value") %>%
          mutate(filtered = input - output) %>%
          pivot_longer(
            cols = c("input", "output", "filtered"),
            names_to = "type", values_to = "value"
          ) %>%
          mutate(
            stage = factor(stage,
              levels = c("amplicon", "demux", "collapse", "graph")
            ),
            M_reads = value / 1e6
          ) %>%
          select(-name)
      }) %>%
      bind_rows(.id = "sample_alias") %>%
      bind_rows(filter(., stage == "amplicon", type == "input") %>%
        mutate(
          stage = "fastQ",
          type = "output"
        )) %>%
      mutate(stage = factor(stage,
        levels = c(
          "fastQ", "amplicon", "demux",
          "collapse", "graph"
        )
      )) %>%
      filter(type != "input") %>%
      arrange(stage, type) %>%
      select(-value) %>%
      group_by(sample_alias)

    if (!is.null(sample_levels)) {
      plot_data <-
        plot_data %>%
        mutate(sample_alias = factor(sample_alias, sample_levels))
    }

    plot_range <-
      max(plot_data$M_reads)

    plots <-
      plot_data %>%
      group_split() %>%
      set_names(group_keys(plot_data)$sample_alias) %>%
      lapply(function(sample_data) {
        sample_data %>%
          ggplot(aes(stage, M_reads, fill = type)) +
          geom_col(show.legend = FALSE) +
          geom_text(
            data = . %>%
              mutate(label = case_when(
                type == "filtered" & stage == "amplicon" ~
                  paste0("Filtered: \n", round(M_reads, 2), " M reads"),
                type == "filtered" & stage == "demux" ~
                  paste0("Demuxed: \n", round(M_reads, 2), " M reads"),
                type == "filtered" & stage == "collapse" ~
                  paste0("Collapsed: \n", round(M_reads, 2), " M reads"),
                type == "filtered" & stage == "graph" ~
                  paste0("Filtered: \n", round(M_reads, 2), " M UEIs"),
                TRUE ~ ""
              )),
            aes(unclass(factor(stage)) - 0.45, label = label),
            position = position_stack(vjust = 0),
            size = 3,
            hjust = 0,
            vjust = 1.1
          ) +
          facet_grid(sample_alias ~ .) +
          scale_y_continuous(
            expand = expansion(c(0, 0.2)),
            limits = c(0, plot_range)
          ) +
          scale_fill_manual(values = c(
            "output" = "#DAD6D7",
            "filtered" = "#C86584"
          )) +
          theme_bw() +
          theme(
            axis.text.x = element_text(angle = 60, hjust = 1),
            panel.border = element_rect(color = "#DAD6D7"),
            panel.grid = element_blank()
          ) +
          labs(x = NULL, y = "Million reads/molecules")
      })


    tabl <-
      plot_data %>%
      mutate(M_reads = round(M_reads, 3)) %>%
      pivot_wider(names_from = c("type"), values_from = "M_reads") %>%
      arrange(sample_alias) %>%
      mutate(`% filtered` = paste0(round(filtered / (output + filtered) * 100, 2), " %")) %>%
      style_table(caption = "Reads/Proteins removed per stage", interactive = FALSE)

    return(list(
      plots = plots,
      table = tabl
    ))
  }

#' Create the component reads per cell
#'
#' This function creates a component that visualizes the reads per cell per sample.
#'
#' @param object A Seurat object.
#'
#' @return A list containing a plot and a table.
#'
#' @export
#'
component_sequencing_reads_per_cell <-
  function(object) {
    plot_data <-
      FetchData(object, c("sample_alias", "reads_in_component"))

    p <-
      plot_data %>%
      plot_violin(
        x = "sample_alias",
        y = "reads_in_component",
        y_label = "Reads per cell",
        round = 1,
        use_1k = TRUE
      )

    tabl <-
      plot_data %>%
      group_by(sample_alias) %>%
      summarise(median_reads_in_component = median(reads_in_component, na.rm = TRUE)) %>%
      select(
        `Sample ID` = sample_alias,
        `Median reads per cell [thousands]` = median_reads_in_component
      ) %>%
      mutate(`Median reads per cell [thousands]` = round(`Median reads per cell [thousands]` / 1e3, 3)) %>%
      style_table(caption = "Median reads per cell", interactive = FALSE)

    return(list(
      plot = p,
      table = tabl
    ))
  }

#' Create the component cell recovery
#'
#' This function creates a component that visualizes the cell recovery
#' after filtering components based on size.
#'
#' @param sample_qc_metrics A list of sample QC metrics.
#' @param sample_levels Optional vector of sample levels to order the samples in the plots.
#'
#' @return A list containing a plot and a table summarizing the cell recovery.
#'
#' @export
#'
component_cell_recovery <-
  function(sample_qc_metrics, sample_levels = NULL) {
    plot_data1 <-
      extract_sample_qc_metrics(
        sample_qc_metrics,
        c(
          component_n_pre_filtering =
            "component_count_pre_component_size_filtering",
          component_n_post_filtering =
            "component_count_post_component_size_filtering",
          min_size_theshold =
            "component_size_min_filtering_threshold",
          max_size_theshold =
            "component_size_max_filtering_threshold"
        ),
        "graph"
      )


    if (!is.null(sample_levels)) {
      plot_data1 <-
        plot_data1 %>%
        mutate(sample_alias = factor(sample_alias, sample_levels))
    }

    p1 <-
      plot_data1 %>%
      select(
        sample_alias,
        component_n_pre_filtering,
        component_n_post_filtering
      ) %>%
      pivot_longer(
        cols = c("component_n_pre_filtering", "component_n_post_filtering"),
        names_to = "type", values_to = "n"
      ) %>%
      mutate(type = str_extract(type, "pre|post") %>%
        str_to_sentence() %>%
        paste("filtering") %>%
        factor(c("Pre filtering", "Post filtering"))) %>%
      ggplot(aes(sample_alias, n)) +
      geom_col(fill = "#DAD6D7") +
      geom_text(aes(label = n),
        vjust = -.1,
        size = 3
      ) +
      facet_wrap(~type, scales = "free", nrow = 1) +
      scale_y_continuous(expand = expansion(c(0, 0.2))) +
      theme_bw() +
      theme(
        axis.text.x = element_text(angle = 60, hjust = 1),
        panel.grid = element_blank()
      ) +
      labs(x = NULL, y = "Number of components")


    plot_data2 <-
      sample_qc_metrics %>%
      map(. %>%
        {
          .$graph$pre_filtering_component_sizes
        } %>%
        unlist() %>%
        enframe("nodes", "n") %>%
        mutate(nodes = as.integer(nodes))) %>%
      bind_rows(.id = "sample_alias") %>%
      filter(nodes >= 20) %>%
      group_by_all() %>%
      reframe(rep = seq_len(n)) %>%
      arrange(sample_alias, -nodes) %>%
      group_by(sample_alias) %>%
      mutate(rank = row_number()) %>%
      ungroup()


    if (!is.null(sample_levels)) {
      plot_data2 <-
        plot_data2 %>%
        mutate(sample_alias = factor(sample_alias, sample_levels))
    }


    p2 <-
      plot_data2 %>%
      ggplot(aes(rank, nodes)) +
      geom_point(size = 0.1) +
      geom_hline(
        data = plot_data1,
        aes(yintercept = min_size_theshold),
        linetype = "dashed"
      ) +
      geom_text(
        data = plot_data1,
        aes(
          x = 1,
          y = min_size_theshold,
          label = min_size_theshold
        ),
        vjust = -0.5,
        hjust = 0
      ) +
      geom_hline(
        data = plot_data1,
        aes(yintercept = max_size_theshold),
        linetype = "dashed"
      ) +
      facet_wrap(~sample_alias) +
      scale_x_log10() +
      scale_y_log10() +
      theme_bw() +
      labs(x = "Component rank", y = "Component size (#nodes)")

    tabl1 <-
      plot_data1 %>%
      arrange(sample_alias) %>%
      select(
        `Sample ID` = sample_alias,
        `Min size threshold` = min_size_theshold,
        `Max size threshold` = max_size_theshold
      ) %>%
      style_table(caption = "Component size thresholds", interactive = FALSE)

    tabl2 <-
      plot_data1 %>%
      arrange(sample_alias) %>%
      select(
        `Sample ID` = sample_alias,
        `Before filtering` = component_n_pre_filtering,
        `After filtering` = component_n_post_filtering
      ) %>%
      style_table(caption = "Number of components", interactive = FALSE)

    return(list(
      plot = list(p1, p2),
      table = list(tabl1, tabl2)
    ))
  }


#' Create the component node and edge count
#'
#' This function creates a component that visualizes the number of nodes and edges
#' before and after cell recovery.
#'
#' @param sample_qc_metrics A list of sample QC metrics.
#' @param sample_levels Optional vector of sample levels to order the samples in the plots.
#'
#' @return A list containing two plots and a table summarizing the node and edge counts.
#'
#' @export
#'
component_node_edge_count <-
  function(sample_qc_metrics, sample_levels = NULL) {
    plot_data <-
      extract_sample_qc_metrics(
        sample_qc_metrics,
        c(
          node_count_pre_recovery = "node_count_pre_recovery",
          node_count_post_recovery = "node_count_post_recovery",
          edge_count_pre_recovery = "edge_count_pre_recovery",
          edge_count_post_recovery = "edge_count_post_recovery"
        ),
        "graph"
      )

    if (!is.null(sample_levels)) {
      plot_data <-
        plot_data %>%
        mutate(sample_alias = factor(sample_alias, sample_levels))
    }


    p1 <-
      plot_data %>%
      select(sample_alias, node_count_pre_recovery, node_count_post_recovery) %>%
      pivot_longer(
        cols = c("node_count_pre_recovery", "node_count_post_recovery"),
        names_to = "type", values_to = "n"
      ) %>%
      mutate(type = str_remove(type, "node_count_") %>%
        str_replace("_", " ") %>%
        str_to_sentence() %>%
        factor(c("Pre recovery", "Post recovery"))) %>%
      ggplot(aes(type, n, fill = type)) +
      geom_col(position = position_dodge()) +
      geom_text(aes(label = paste(round(n / 1e6, 3), " M")),
        vjust = -.1,
        size = 3
      ) +
      facet_wrap(~sample_alias, nrow = 1) +
      scale_y_continuous(expand = expansion(c(0, 0.2))) +
      scale_fill_manual(values = c("#DAD6D7", "#B4ADAF")) +
      theme_bw() +
      theme(
        axis.text.x = element_text(angle = 60, hjust = 1),
        panel.grid = element_blank()
      ) +
      labs(x = NULL, y = "Number of nodes")

    p2 <-
      plot_data %>%
      select(sample_alias, edge_count_pre_recovery, edge_count_post_recovery) %>%
      pivot_longer(
        cols = c("edge_count_pre_recovery", "edge_count_post_recovery"),
        names_to = "type", values_to = "n"
      ) %>%
      mutate(type = str_remove(type, "edge_count_") %>%
        str_replace("_", " ") %>%
        str_to_sentence() %>%
        factor(c("Pre recovery", "Post recovery"))) %>%
      ggplot(aes(type, n, fill = type)) +
      geom_col(position = position_dodge()) +
      geom_text(aes(label = paste(round(n / 1e6, 3), " M")),
        vjust = -.1,
        size = 3
      ) +
      facet_wrap(~sample_alias, nrow = 1) +
      scale_y_continuous(expand = expansion(c(0, 0.2))) +
      scale_fill_manual(values = c("#DAD6D7", "#B4ADAF")) +
      theme_bw() +
      theme(
        axis.text.x = element_text(angle = 60, hjust = 1),
        panel.grid = element_blank()
      ) +
      labs(x = NULL, y = "Number of edges")

    tabl <-
      plot_data %>%
      arrange(sample_alias) %>%
      select(
        `Sample ID` = sample_alias,
        `Node count before recovery` = node_count_pre_recovery,
        `Node count after recovery` = node_count_post_recovery,
        `Edge count before recovery` = edge_count_pre_recovery,
        `Edge count after recovery` = edge_count_post_recovery
      ) %>%
      mutate(
        `Node count before recovery` = round(`Node count before recovery` /
          1e6, 3),
        `Node count after recovery` = round(`Node count after recovery` /
          1e6, 3),
        `Edge count before recovery` = round(`Edge count before recovery` /
          1e6, 3),
        `Edge count after recovery` = round(`Edge count after recovery` /
          1e6, 3)
      ) %>%
      style_table(caption = "Node and edge counts [millions]", interactive = FALSE)

    return(list(
      plots = list(p1, p2),
      table = tabl
    ))
  }

#' Create the component node degree
#'
#' This function creates a component that visualizes the mean node degree for A and B nodes
#' for each sample.
#'
#' @param object A Seurat object containing the sample data.
#' @param sample_levels Optional vector of sample levels to order the samples in the plots.
#'
#' @return A list containing a plot and a table summarizing the mean node degree.
#'
#' @export
#'
component_node_degree <-
  function(object, sample_levels = NULL) {
    plot_data <-
      FetchData(
        object,
        c(
          "sample_alias",
          "A_nodes_mean_degree",
          "B_nodes_mean_degree"
        )
      ) %>%
      as_tibble(rownames = "sample_component") %>%
      pivot_longer(
        cols = c("A_nodes_mean_degree", "B_nodes_mean_degree"),
        names_to = "type", values_to = "degree"
      ) %>%
      mutate(type = str_remove(type, "_nodes_mean_degree"))


    if (!is.null(sample_levels)) {
      plot_data <-
        plot_data %>%
        mutate(sample_alias = factor(sample_alias, sample_levels))
    }

    p <-
      plot_data %>%
      plot_violin(
        x = "sample_alias",
        y = "degree",
        fill = "type",
        title = "Mean A and B node degree",
        y_label = "Mean node degree",
        fill_label = "Node type",
        palette = yes_no_palette
      )

    tabl <-
      plot_data %>%
      group_by(sample_alias, type) %>%
      summarise(median_degree = median(degree), .groups = "drop") %>%
      pivot_wider(names_from = "type", values_from = "median_degree") %>%
      select(
        `Sample ID` = sample_alias,
        `COa` = A,
        `COb` = B
      ) %>%
      mutate(
        `COa` = round(`COa`, 2),
        `COb` = round(`COb`, 2)
      ) %>%
      style_table(caption = "Median A and B node degree", interactive = FALSE)

    return(list(
      plot = p,
      table = tabl
    ))
  }

#' Create the component crossing edges
#'
#' This function creates a component that visualizes the percentage of crossing edges
#' for each sample, both in the initial and refinement stages.
#'
#' @param qc_metrics_tables A list of QC metrics tables.
#' @param sample_levels Optional vector of sample levels to order the samples in the plots.
#'
#' @return A list containing a plot and a table summarizing the crossing edges.
#'
#' @export
#'
component_crossing_edges <-
  function(qc_metrics_tables, sample_levels = NULL) {
    plot_data <-
      qc_metrics_tables$crossing_edges %>%
      mutate(label = paste0(edges, "\n", round(percent, 3), " %"))

    if (!is.null(sample_levels)) {
      plot_data <-
        plot_data %>%
        mutate(sample_alias = factor(sample_alias, sample_levels))
    }
    p <-
      plot_data %>%
      ggplot(aes(sample_alias, percent, fill = type)) +
      geom_col() +
      geom_text(
        data = . %>%
          mutate(label = ifelse(type == "Initial stage",
            label, ""
          )),
        aes(
          x = unclass(factor(sample_alias)) - 0.45,
          label = label
        ),
        position = position_stack(vjust = 1),
        vjust = 1,
        hjust = 0,
        size = 3
      ) +
      geom_text(
        data = . %>%
          mutate(label = ifelse(type == "Refinement stage",
            label, ""
          )),
        aes(
          x = unclass(factor(sample_alias)) - 0.45,
          label = label
        ),
        position = position_stack(vjust = 1),
        vjust = -0.1,
        hjust = 0,
        size = 3
      ) +
      geom_text(
        aes(
          x = unclass(factor(sample_alias)) + 0.45,
          y = 100 * removed_total / total_edges_in,
          label = paste0(
            removed_total, "\n",
            round(100 * removed_total / total_edges_in, 3), " %"
          )
        ),
        vjust = -0.1,
        hjust = 1,
        fontface = "bold",
        size = 3
      ) +
      scale_y_continuous(expand = expansion(c(0, 0.2))) +
      scale_fill_manual(values = c("#DAD6D7", "#B4ADAF")) +
      theme_bw() +
      theme(
        axis.text.x = element_text(angle = 60, hjust = 1),
        panel.grid = element_blank()
      ) +
      labs(x = NULL, y = "Crossing edges removed [%]")

    tabl <-
      qc_metrics_tables$crossing_edges %>%
      arrange(sample_alias) %>%
      mutate(
        percent = paste(round(percent, 3), "%"),
        percent_total = paste(
          round(100 * removed_total / total_edges_in, 3),
          "%"
        )
      ) %>%
      pivot_wider(names_from = "type", values_from = c("edges", "percent")) %>%
      unite("Initial stage",
        c("edges_Initial stage", "percent_Initial stage"),
        sep = " / "
      ) %>%
      unite("Refinement stage",
        c("edges_Refinement stage", "percent_Refinement stage"),
        sep = " / "
      ) %>%
      unite("Total removed",
        c("removed_total", "percent_total"),
        sep = " / "
      ) %>%
      select(-total_edges_in) %>%
      style_table(caption = "Crossing edges removed", interactive = FALSE)

    return(list(
      plot = p,
      table = tabl
    ))
  }

#' Create the component coreness
#'
#' This function creates a component that visualizes the k-coreness of nodes in each component.
#'
#' @param qc_metrics_tables A list of QC metrics tables containing coreness data.
#'
#' @return A list containing three plots and a table summarizing the coreness metrics.
#'
#' @export
#'
component_coreness <-
  function(qc_metrics_tables) {
    p1 <-
      qc_metrics_tables$coreness$component_summary %>%
      plot_violin(
        x = "sample_alias",
        y = "mean_coreness",
        title = "Mean coreness",
        subtitle = "The average coreness of nodes in each component",
        y_label = "Mean coreness"
      )

    # Dangling nodes - % coreness 1
    p2 <-
      qc_metrics_tables$coreness$component_summary %>%
      plot_violin(
        x = "sample_alias",
        y = "percent_dangling_nodes",
        title = "Percent dangling nodes",
        subtitle = "The percentage of nodes with coreness 1",
        y_label = "Dangling nodes [%]", round = 2
      )

    # Well-connected nodes - % coreness 3-5
    p3 <-
      qc_metrics_tables$coreness$component_summary %>%
      plot_violin(
        x = "sample_alias",
        y = "percent_well_connected_nodes",
        title = "Percent well-connected nodes",
        subtitle = "The percentage of nodes with coreness 3-5",
        y_label = "Well-connected nodes [%]",
        round = 2
      )

    tabl <-
      qc_metrics_tables$coreness$sample_summary %>%
      select(
        `Sample ID` = sample_alias,
        `Median mean coreness` = median_mean_coreness,
        `Median dangling nodes [%]` = median_percent_dangling_nodes,
        `Median well connected nodes [%]` = median_percent_well_connected_nodes
      ) %>%
      mutate(
        `Median mean coreness` =
          round(`Median mean coreness`, 2),
        `Median dangling nodes [%]` =
          round(`Median dangling nodes [%]`, 2),
        `Median well connected nodes [%]` =
          round(`Median well connected nodes [%]`, 2)
      ) %>%
      style_table(caption = "k-coreness", interactive = FALSE)

    heatmap_data <-
      qc_metrics_tables$coreness$data %>%
      group_by(sample_alias)

    heatmap_plots <-
      heatmap_data %>%
      group_split() %>%
      set_names(group_keys(heatmap_data)$sample_alias) %>%
      lapply(function(g_data) {
        if (n_distinct(g_data$sample_component) > 2000) {
          plot_cells <-
            g_data %>%
            select(sample_component, n_umi) %>%
            distinct() %>%
            arrange(-n_umi, sample_component) %>%
            pull(sample_component) %>%
            head(2000)

          g_data <-
            g_data %>%
            filter(sample_component %in% plot_cells)
        }

        g_data %>%
          select(sample_component, k_core, percent_nodes) %>%
          pivot_wider(names_from = "k_core", values_from = "percent_nodes") %>%
          column_to_rownames("sample_component") %>%
          as.matrix() %>%
          t() %>%
          ComplexHeatmap::pheatmap(
            color = cherry_gradient,
            breaks = seq(0, 100, 1),
            cellheight = 10,
            show_colnames = FALSE,
            cluster_rows = FALSE,
            clustering_method = "ward.D2",
            heatmap_legend_param = list(title = "% Nodes")
          ) %>%
          as.ggplot()
      })

    return(list(
      plots = list(p1, p2, p3),
      table = tabl,
      heatmap_plots = heatmap_plots
    ))
  }

#' Create the component for abundance per sample
#'
#' This function creates plots visualizing the abundance of each marker
#' across different samples.
#'
#' @param object A processed data object containing marker abundance data.
#' @param params A list of parameters, including control markers.
#' @param metadata A metadata object containing sample information.
#' @param sample_palette A color palette for the samples.
#' @param test_mode A boolean indicating whether to run in test mode (default is FALSE).
#'
#' @return A list containing plots.
#'
#' @export
#'
component_abundance_per_sample <- function(
  object,
  params,
  metadata,
  sample_palette,
  test_mode = FALSE
) {
  plot_data <-
    object %>%
    LayerData("data") %>%
    as_tibble(rownames = "marker") %>%
    {
      if (test_mode) {
        slice_head(., n = 10)
      } else {
        .
      }
    } %>%
    pivot_longer(-marker, names_to = "cell_id", values_to = "normcount") %>%
    left_join(
      FetchData(object,
        vars = c("sample_alias", "condition")
      ) %>%
        as_tibble(rownames = "cell_id"),
      by = "cell_id"
    ) %>%
    mutate(
      marker = factor(marker, order_cd_markers(
        unique(marker),
        params$control_markers
      )),
      sample_alias = factor(sample_alias, levels = metadata$sample_alias)
    ) %>%
    group_by(marker)

  plots <-
    plot_data %>%
    group_split() %>%
    set_names(group_keys(plot_data)$marker) %>%
    lapply(function(g_data) {
      g_data %>%
        plot_violin(
          x = "sample_alias",
          y = "normcount",
          fill = "condition",
          fill_label = "Condition",
          title = unique(g_data$marker) %>% as.character(),
          y_label = "Normalized counts",
          summarize = FALSE,
          palette = sample_palette,
          use_jitter = TRUE,
          use_grid = TRUE,
          jitter_alpha = 1
        )
    })

  return(plots)
}


#' Create the component for abundance per celltype
#'
#' This function creates plots visualizing the abundance of each marker
#' across different celltypes
#'
#' @param object A processed data object containing marker abundance data.
#' @param params A list of parameters containing control markers.
#' @param sample_palette A color palette for the samples.
#' @param test_mode A boolean indicating whether to run in test mode (default is FALSE).
#'
#' @return A list containing plots.
#'
#' @export
#'
component_abundance_per_celltype <- function(
  object,
  params,
  sample_palette,
  test_mode = FALSE
) {
  plot_data <-
    object %>%
    LayerData("data") %>%
    as_tibble(rownames = "marker") %>%
    {
      if (test_mode) {
        slice_head(., n = 10)
      } else {
        .
      }
    } %>%
    pivot_longer(-marker, names_to = "cell_id", values_to = "normcount") %>%
    left_join(
      FetchData(object,
        vars = c(
          "sample_alias",
          "seurat_clusters",
          "l1_annotation_summary"
        )
      ) %>%
        as_tibble(rownames = "cell_id"),
      by = "cell_id"
    ) %>%
    mutate(
      marker = factor(marker, order_cd_markers(
        unique(marker),
        params$control_markers
      )),
      l1_annotation_summary = factor(l1_annotation_summary)
    ) %>%
    group_by(marker)

  plots <-
    plot_data %>%
    group_split() %>%
    set_names(group_keys(plot_data)$marker) %>%
    lapply(function(g_data) {
      g_data %>%
        plot_violin(
          x = "sample_alias",
          y = "normcount",
          fill = "sample_alias",
          title = unique(g_data$marker) %>% as.character(),
          y_label = "Normalized counts",
          summarize = FALSE,
          palette = sample_palette,
          facet_var = "l1_annotation_summary",
          use_jitter = TRUE,
          use_grid = TRUE,
          jitter_alpha = 1
        ) +
        guides(fill = "none")
    })

  return(plots)
}


#' Create the component for proximity Z scores
#'
#' This function creates plots visualizing the proximity Z scores for selected
#' contrasts.
#'
#' @param processed_data A processed data object containing marker abundance data.
#' @param sample_palette A color palette for the samples.
#' @param proximity_score One of "join_count_z" or "log2_ratio".
#' @param test_mode A boolean indicating whether to run in test mode (default is FALSE).
#'
#' @return A list containing plots.
#'
#' @export
#'
component_proximity_selected <- function(
  processed_data,
  sample_palette,
  proximity_score,
  test_mode = FALSE
) {
  plots <-
    processed_data %>%
    group_split() %>%
    set_names(group_keys(processed_data)$contrast) %>%
    {
      if (test_mode) {
        head(., n = 3)
      } else {
        .
      }
    } %>%
    lapply(function(g_data) {
      g_data %>%
        complete(
          sample_alias = levels(g_data$sample_alias),
          l1_annotation_summary = c("CD4 T", "CD8 T", "NK", "Mono", "B"),
          fill = setNames(list(NA), proximity_score)
        ) %>%
        ggplot(aes(x = sample_alias, y = !!sym(proximity_score))) +
        geom_hline(yintercept = 0) +
        geom_violin(aes(fill = condition), draw_quantiles = 0.5, drop = FALSE) +
        geom_jitter(
          size = 0.1,
          position = position_jitter(height = 0),
          alpha = 0.3
        ) +
        geom_text(
          data = . %>%
            group_by(sample_alias, l1_annotation_summary) %>%
            summarise(!!sym(proximity_score) := median(!!sym(proximity_score)), .groups = "drop"),
          aes(label = round(!!sym(proximity_score), 2)),
          size = 3,
          vjust = 0
        ) +
        facet_grid(l1_annotation_summary ~ .) +
        scale_fill_manual(values = sample_palette) +
        theme_violin() +
        theme(
          axis.text.x = element_text(angle = 45, hjust = 1),
          panel.grid = element_blank()
        ) +
        labs(
          title = g_data$contrast[1],
          x = "Sample ID",
          y = ifelse(proximity_score == "join_count_z",
            "Z-score",
            "logratio"
          ),
        )
    })

  return(plots)
}


#' Create the component for proximity per marker
#'
#' This function creates plots visualizing the proximity scores for each marker
#' across different samples and conditions.
#'
#' @param processed_data A processed data object containing proximity scores.
#' @param proximity_score One of "join_count_z" or "log2_ratio".
#' @param sample_palette A color palette for the samples.
#' @param test_mode A boolean indicating whether to run in test mode (default is FALSE).
#'
#' @return A list containing plots for each marker.
#'
#' @export
#'
component_proximity_per_marker <- function(
  processed_data,
  proximity_score,
  sample_palette,
  test_mode = FALSE) {
  plots <-
    processed_data %>%
    group_split() %>%
    set_names(group_keys(processed_data)$marker_1) %>%
    {
      if (test_mode) {
        head(., 10)
      } else {
        .
      }
    } %>%
    lapply(function(g_data) {
      g_data %>%
        complete(
          sample_alias = levels(g_data$sample_alias),
          condition = levels(g_data$condition),
          fill = setNames(list(NA), proximity_score)
        ) %>%
        ggplot(aes(sample_alias, !!sym(proximity_score), fill = condition)) +
        geom_hline(yintercept = 0) +
        geom_violin(
          draw_quantiles = 0.5,
          scale = "width"
        ) +
        geom_jitter(
          size = 0.1,
          position = position_jitter(height = 0)
        ) +
        scale_fill_manual(values = sample_palette) +
        theme_violin() +
        labs(
          x = "Sample ID",
          y = ifelse(
            proximity_score == "join_count_z",
            "Self-Proximity Z-score",
            "Self-Proximity logratio"
          ),
          title = unique(g_data$marker_1)
        )
    })
}

#' Create a clustering summary component
#'
#' This function generates summary plots of the clustering scores.
#'
#' @param proximity_scores A data frame containing proximity scores.
#' @param heatmap_gradient A color gradient for the heatmap.
#' @param test_mode A boolean indicating whether to run in test mode (default is FALSE).
#'
#' @return A list of ggplot objects summarizing the clustering scores.
#'
#' @export
#'
component_clustering_summary <- function(
  proximity_scores,
  heatmap_gradient,
  test_mode = FALSE
) {
  plot_data <-
    proximity_scores %>%
    filter(marker_1 == marker_2) %>%
    rename(component = sample_component) %>%
    SummarizeProximityScores(
      proximity_metric = "join_count_z", detailed = TRUE,
      group_vars = c(
        "sample_alias",
        "l1_annotation_summary",
        "condition"
      )
    ) %>%
    rowwise() %>%
    mutate(mean_log2_ratio = mean(log2(
      pmax(unlist(join_count_list), 1) /
        pmax(unlist(join_count_expected_mean_list), 1)
    ))) %>%
    rename(marker = marker_1) %>%
    {
      if (test_mode) {
        .
      } else {
        filter(., n_cells_detected >= 10)
      }
    } %>%
    group_by(marker) %>%
    filter(any(abs(mean_join_count_z) > 1)) %>%
    ungroup() %>%
    select(
      sample_alias, condition, marker, l1_annotation_summary,
      mean_join_count_z, mean_log2_ratio
    )

  plot_join_count_z <-
    plot_data %>%
    ggplot(aes(marker, l1_annotation_summary, fill = mean_join_count_z)) +
    geom_tile() +
    facet_grid(condition ~ ., scales = "free", space = "free") +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
      legend.position = "bottom"
    ) +
    scale_fill_gradientn(
      colors = heatmap_gradient,
      limits = c(-1, 1) * max(abs(plot_data$mean_join_count_z))
    ) +
    labs(
      fill = "Mean Z-score",
      x = "",
      y = "Celltype"
    )

  plot_log2_ratio <-
    plot_data %>%
    ggplot(aes(marker, l1_annotation_summary, fill = mean_log2_ratio)) +
    geom_tile() +
    facet_grid(condition ~ ., scales = "free", space = "free") +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
      legend.position = "bottom"
    ) +
    scale_fill_gradientn(
      colors = heatmap_gradient,
      limits = c(-1, 1) * max(abs(plot_data$mean_log2_ratio))
    ) +
    labs(
      fill = "Mean logratio",
      x = "",
      y = "Celltype"
    )

  return(list("join_count_z" = plot_join_count_z, "log2_ratio" = plot_log2_ratio))
}


#' Create the component for proximity heatmaps per sample and condition
#'
#' This function creates plots visualizing the proximity scores as heatmaps per sample and condition.
#'
#' @param processed_data A processed data object containing marker abundance data.
#' @param heatmap_gradient A color palette for the heatmaps.
#' @param proximity_score One of "join_count_z" or "log2_ratio".
#'
#' @return A list containing plots.
#'
#' @export
#'
component_proximity_heatmap_sample <- function(
  processed_data,
  heatmap_gradient,
  proximity_score
) {
  plots <-
    processed_data %>%
    group_split() %>%
    set_names(group_keys(processed_data)$sample_alias) %>%
    lapply(function(g_data) {
      sampl <- unique(g_data$sample_alias)
      cond <- unique(g_data$condition)

      plot_markers <-
        g_data %>%
        bind_rows(select(g_data,
          marker_2 = marker_1,
          marker_1 = marker_2,
          everything()
        )) %>%
        distinct() %>%
        group_by(marker = marker_1) %>%
        summarise(sd = sd(mean_join_count_z)) %>%
        arrange(-sd) %>%
        head(60) %>%
        pull(marker)

      g_data %>%
        filter(marker_1 %in% plot_markers & marker_2 %in% plot_markers) %>%
        ColocalizationHeatmap(
          marker1_col = "marker_1",
          marker2_col = "marker_2",
          value_col = ifelse(
            proximity_score == "join_count_z",
            "mean_join_count_z",
            "mean_log2_ratio"
          ),
          clustering_method = "ward.D2",
          cellwidth = 8,
          cellheight = 8,
          fontsize = 8,
          colors = heatmap_gradient,
          border_color = "white",
          legend_range = if (proximity_score == "join_count_z") {
            c(-1, 1) * 4
          } else {
            c(-1, 1)
          },
          legend_title = ifelse(
            proximity_score == "join_count_z",
            "Mean\nProximity\nZ-score",
            "Mean\nProximity\nlogratio"
          ),
          main = paste("Mean colocalization:", sampl, "-", cond)
        ) %>%
        as.ggplot()
    })

  return(plots)
}


#' Create the component for proximity heatmaps per celltype and condition
#'
#' This function creates plots visualizing the proximity scores as heatmaps per celltype and condition.
#'
#' @param processed_data A processed data object containing marker abundance data.
#' @param heatmap_gradient A color palette for the heatmaps.
#' @param proximity_score One of "join_count_z" or "log2_ratio".
#'
#' @return A list containing plots.
#'
#' @export
#'
component_proximity_heatmap_celltype <- function(
  processed_data,
  heatmap_gradient,
  proximity_score
) {
  plots <-
    processed_data %>%
    group_split() %>%
    set_names(paste(
      group_keys(processed_data)$celltype,
      group_keys(processed_data)$condition
    )) %>%
    lapply(function(g_data) {
      cellpop <- unique(g_data$celltype)
      cond_ <- unique(g_data$condition)

      plot_markers <-
        g_data %>%
        bind_rows(select(g_data,
          marker_2 = marker_1,
          marker_1 = marker_2,
          everything()
        )) %>%
        distinct() %>%
        group_by(marker = marker_1) %>%
        summarise(sd = sd(mean_join_count_z)) %>%
        arrange(-sd) %>%
        head(60) %>%
        pull(marker)


      g_data %>%
        filter(marker_1 %in% plot_markers & marker_2 %in% plot_markers) %>%
        ColocalizationHeatmap(
          marker1_col = "marker_1",
          marker2_col = "marker_2",
          value_col = ifelse(
            proximity_score == "join_count_z",
            "mean_join_count_z",
            "mean_log2_ratio"
          ),
          clustering_method = "ward.D2",
          cellwidth = 8,
          cellheight = 8,
          fontsize = 8,
          colors = heatmap_gradient,
          border_color = "white",
          legend_range = if (proximity_score == "join_count_z") {
            c(-1, 1) * 4
          } else {
            c(-1, 1)
          },
          legend_title = ifelse(
            proximity_score == "join_count_z",
            "Mean\nProximity\nZ-score",
            "Mean\nProximity\nlogratio"
          ),
          main = paste("Mean colocalization:", cond_, cellpop)
        ) %>%
        as.ggplot()
    })

  return(plots)
}

#' Create the component dimensionality reduction plots
#'
#' This function creates a component that visualizes the dimensionality reduction plots
#' for PCA, colored by clusters, samples, and conditions.
#'
#' @param object A Seurat object containing the processed data.
#' @param plot_reductions A character vector specifying which dimensionality reductions to plot.
#' @param sample_palette Optional color palette for samples.
#' @param cluster_palette Optional color palette for clusters.
#' @param sample_plots Logical indicating whether to include sample-wise plots.
#'
#' @return A list containing three plots for reductions.
#'
#' @export
component_dimred_plots <-
  function(
    object,
    plot_reductions = "pca",
    sample_palette = NULL,
    cluster_palette = NULL,
    sample_plots = FALSE
  ) {
    plot_reductions %>%
      set_names(toupper(plot_reductions)) %>%
      lapply(function(red) {
        p1 <-
          plot_embedding(object, red,
            metavars = "seurat_clusters",
            pal = cluster_palette,
            plot_title = "Clusters",
            xaxis_title = FALSE
          )
        p2 <-
          plot_embedding(object, red,
            metavars = "sample_alias",
            pal = sample_palette,
            plot_title = "Samples",
            label = FALSE,
            legend_position = "bottom",
            yaxis_title = FALSE
          )
        p3 <-
          plot_embedding(object, red,
            metavars = "condition",
            pal = sample_palette,
            plot_title = "Conditions",
            label = FALSE,
            legend_position = "bottom",
            xaxis_title = FALSE,
            yaxis_title = FALSE
          )

        plot_list <-
          list(Combined = p1 | p2 | p3)

        if (sample_plots) {
          samplewise_plots <-
            plot_embeddings_samplewise(
              levels(object$sample_alias),
              object,
              plot_reduction = red,
              pal = cluster_palette
            )

          plot_list <-
            list(
              Combined = p1 | p2 | p3,
              Samplewise = samplewise_plots
            )
        } else {
          plot_list <-
            p1 | p2 | p3
        }

        return(plot_list)
      })
  }


#' Create the component annotation summary
#'
#' This function creates a component that visualizes the annotation summary
#' for l1 and l2 annotations in the specified reduction.
#'
#' @param object A Seurat object containing the processed data.
#' @param heatmap_gradient A color gradient for the heatmap.
#' @param reduction Optional character specifying the reduction to use. If NULL, the first preferred reduction is used.
#'
#' @return A list of ggplot objects visualizing the l1 and l2 annotation summaries.
#'
#' @export
#'
component_annotation <-
  function(
    object,
    heatmap_gradient,
    reduction = NULL
  ) {
    if (is.null(reduction)) {
      reduction <- preferred_dimred_order(Reductions(object))[1]
    }
    p1 <-
      plot_embedding(
        object,
        reduction,
        metavars = "l1_annotation_summary",
        pal = cell_palette,
        label = TRUE,
        legend_position = "none",
        xaxis_title = TRUE,
        yaxis_title = TRUE
      )
    p2 <-
      plot_embedding(
        object,
        reduction,
        metavars = "l2_annotation_summary",
        pal = cell_palette,
        label = TRUE,
        legend_position = "none",
        xaxis_title = TRUE,
        yaxis_title = TRUE
      )


    plot_data <-
      object %>%
      LayerData("data") %>%
      as_tibble(rownames = "marker") %>%
      pivot_longer(-marker, names_to = "cell_id", values_to = "normcount") %>%
      left_join(FetchData(object,
        vars = c("seurat_clusters", "l1_annotation_summary")
      ) %>%
        as_tibble(rownames = "cell_id")) %>%
      group_by(seurat_clusters, cell_annotation = l1_annotation_summary, marker) %>%
      summarise(median = median(normcount)) %>%
      unite(id, seurat_clusters, cell_annotation) %>%
      pivot_wider(names_from = marker, values_from = median, values_fill = 0) %>%
      column_to_rownames("id") %>%
      as.matrix()

    plot_legend <-
      max(abs(plot_data))

    p3 <-
      plot_data %>%
      t() %>%
      pheatmap(
        color = heatmap_gradient,
        breaks = seq(-0, plot_legend, length.out = 100),
        clustering_method = "ward.D2",
        cellwidth = 7,
        cellheight = 7,
        fontsize = 7,
        heatmap_legend_param = list(title = "Median Norm.\nAbundance"),
        main = "Cluster marker abundance"
      ) %>%
      as.ggplot()

    plot_data <-
      FetchData(object,
        vars = c(
          "seurat_clusters",
          "l1_annotation_summary",
          "condition",
          "sample_alias"
        )
      ) %>%
      as_tibble(rownames = "cell_id") %>%
      group_by(sample_alias, condition, l1_annotation_summary) %>%
      count() %>%
      group_by(sample_alias, condition) %>%
      mutate(frac = n / sum(n)) %>%
      ungroup() %>%
      complete(
        nesting(sample_alias, condition),
        l1_annotation_summary,
        fill = list(n = 0, frac = 0)
      ) %>%
      group_by(sample_alias)

    p4 <-
      plot_data %>%
      ggplot(aes(sample_alias, frac, fill = l1_annotation_summary)) +
      geom_col(position = "stack") +
      geom_text(aes(label = scales::percent(frac)),
        position = position_stack(vjust = 0.5), size = 2
      ) +
      scale_fill_manual(values = cell_palette) +
      theme_bw() +
      scale_y_continuous(
        expand = expansion(0),
        labels = scales::percent
      ) +
      theme(
        axis.text.x = element_text(angle = 60, hjust = 1),
        panel.grid = element_blank()
      ) +
      labs(
        x = "Sample",
        y = "% cells",
        fill = "Cell type",
        title = "Cell type composition"
      )

    barplots1 <-
      plot_data %>%
      group_split() %>%
      set_names(group_keys(plot_data)$sample_alias) %>%
      lapply(function(g_data) {
        g_data %>%
          ggplot(aes(l1_annotation_summary, n, fill = l1_annotation_summary)) +
          geom_col() +
          geom_text(aes(label = paste(n, scales::percent(frac), sep = "\n")),
            position = position_dodge(width = 0.9),
            vjust = -0.25,
            size = 2
          ) +
          scale_fill_manual(values = cell_palette) +
          theme_bw() +
          scale_y_continuous(expand = expansion(c(0, 0.2))) +
          theme(
            axis.text.x = element_text(angle = 60, hjust = 1),
            panel.grid = element_blank()
          ) +
          labs(
            x = "Cell type",
            y = "Number of cells",
            fill = "Cell type",
            title = "Number of cells per annotation",
            subtitle = unique(g_data$sample_alias)
          )
      })

    plot_data <-
      plot_data %>%
      arrange(l1_annotation_summary) %>%
      group_by(l1_annotation_summary)

    barplots2 <-
      plot_data %>%
      group_split() %>%
      set_names(group_keys(plot_data)$l1_annotation_summary) %>%
      lapply(function(g_data) {
        g_data %>%
          ggplot(aes(sample_alias, frac)) +
          geom_col(
            show.legend = FALSE,
            fill = "#DAD6D7"
          ) +
          geom_text(aes(label = scales::percent(frac)),
            position = position_stack(vjust = 1),
            vjust = -0.5,
            size = 2
          ) +
          theme_bw() +
          scale_y_continuous(
            expand = expansion(0), limits = c(0, 1),
            labels = scales::percent
          ) +
          theme(
            axis.text.x = element_text(angle = 60, hjust = 1),
            panel.grid = element_blank()
          ) +
          labs(
            x = "Sample",
            y = "% cells",
            title = paste(
              "Fraction of",
              unique(g_data$l1_annotation_summary),
              "cells per sample"
            )
          )
      })

    tabl <-
      plot_data %>%
      ungroup() %>%
      mutate(frac = scales::percent(frac, accuracy = 0.1)) %>%
      select(
        `Sample ID` = sample_alias,
        `Cell annotation` = l1_annotation_summary,
        frac
      ) %>%
      pivot_wider(names_from = `Cell annotation`, values_from = frac) %>%
      style_table(caption = "Cell type composition [%]", interactive = FALSE)

    return(list(
      dimred_plots = list(
        "Level 1 annotation" = p1,
        "Level 2 annotation" = p2
      ),
      heatmap = p3,
      celltype_composition = p4,
      celltype_composition_barplots1 = barplots1,
      celltype_composition_barplots2 = barplots2,
      celltype_composition_table = tabl
    ))
  }
