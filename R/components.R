#' Create the denoising component
#'
#' This function creates plots and tables visualizing denoising metrics,
#' including the fraction of UMIs removed, UMIs removed by method, and
#' reduction in isotype fraction.
#'
#' @param qc_metrics_tables A list from [get_qc_metrics()] containing `denoising`
#'   and `denoising_detail` slots.
#' @param sample_levels Optional vector of sample levels to order samples in plots.
#'
#' @return A list with named `plots` and `tables` for removed UMIs, denoising
#'   by method, and isotype reduction.
#'
#' @export
#'
component_denoising <- function(
  qc_metrics_tables,
  sample_levels = NULL
) {
  pixelatorR:::assert_class(qc_metrics_tables, "list")

  if (is.null(qc_metrics_tables$denoising)) {
    cli::cli_abort("Denoising data is required but {.code qc_metrics_tables$denoising} is NULL.")
  }

  denoising_data <- qc_metrics_tables$denoising

  if ("pool" %in% names(denoising_data)) {
    denoising_data <-
      denoising_data %>%
      rename(sample_alias = pool)
  }

  denoising_data <- set_sample_levels(denoising_data, sample_levels)

  p_removed_umis <-
    denoising_data %>%
    ggplot(aes(sample_alias, percent_umis_denoised)) +
    geom_col(fill = "#DAD6D7") +
    geom_text(aes(label = paste0(round(percent_umis_denoised, 3), " %")),
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
      title = "Fraction of UMIs removed by denoising",
      caption = expression(~ italic(frac("# Discarded UMIs", "# Total UMIs")))
    )

  tabl_removed_umis <-
    denoising_data %>%
    arrange(sample_alias) %>%
    mutate(
      percent_umis_denoised = round(percent_umis_denoised, 3),
      number_of_umis_removed = round(number_of_umis_removed / 1e6, 3)
    ) %>%
    select(
      `Sample ID` = sample_alias,
      `% Denoised UMIs` = percent_umis_denoised,
      `Total denoised UMIs [M]` = number_of_umis_removed
    ) %>%
    style_table(caption = "UMIs removed by denoising", interactive = FALSE)

  by_method_data <-
    qc_metrics_tables$denoising_detail$by_method %>%
    set_sample_levels(sample_levels) %>%
    mutate(type = str_replace(type, "_and_", "&"))

  method_palette <-
    c(
      "ace" = "#DAD6D7",
      "ace&pls" = "#C86584",
      "pls" = "#B4ADAF"
    )

  # Add missing methods to palette
  if (any(!by_method_data$type %in% names(method_palette))) {
    missing_methods <- setdiff(unique(by_method_data$type), names(method_palette))
    add_palette <- c("#496389", "#1F395F", "#4D988D", "#E05573", "#BF9871", "#918F8F")
    method_palette <- c(method_palette, setNames(add_palette[seq_along(missing_methods)], missing_methods))
  }

  p_by_method <-
    by_method_data %>%
    mutate(type = factor(type, levels = names(method_palette))) %>%
    arrange(sample_alias, desc(type)) %>%
    group_by(sample_alias) %>%
    mutate(
      fraction_denoised = umis_denoised / n_umi,
      x_pos = unclass(factor(type)) %% 2,
      y_pos = cumsum(fraction_denoised) - 0.5 * fraction_denoised
    ) %>%
    ggplot(aes(sample_alias, fraction_denoised, fill = type)) +
    geom_col() +
    geom_text(
      aes(
        label = paste0(round(100 * fraction_denoised, 1), " %"),
        x = unclass(sample_alias) + c(-0.45, 0.45)[x_pos + 1],
        y = y_pos,
        hjust = x_pos
      ),
      vjust = -.1,
      size = 3
    ) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    scale_fill_manual(values = method_palette, na.value = "#DAD6D7") +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 60, hjust = 1),
      panel.grid = element_blank()
    ) +
    labs(
      x = NULL,
      y = "Fraction of UMIs denoised",
      fill = "Method",
      title = "Fraction of UMIs denoised by method"
    )

  tabl_by_method <-
    by_method_data %>%
    mutate(
      fraction_denoised = round(100 * umis_denoised / n_umi, 3)
    ) %>%
    arrange(sample_alias, type) %>%
    select(
      `Sample ID` = sample_alias,
      `Denoising method` = type,
      `UMIs denoised` = umis_denoised,
      `Fraction denoised [%]` = fraction_denoised
    ) %>%
    style_table(caption = "UMIs denoised by method", interactive = FALSE)

  isotype_data <-
    qc_metrics_tables$denoising_detail$isotype_reduction %>%
    set_sample_levels(sample_levels)

  p_isotype_reduction <-
    isotype_data %>%
    plot_violin(
      x = "sample_alias",
      y = "isotype_reduction",
      title = "Reduction in isotype fraction",
      y_label = "Isotype reduction by denoising",
      use_pct = TRUE,
      hline = 0,
      round = 1
    ) +
    scale_y_continuous(
      labels = scales::percent_format(accuracy = 1),
      limits = c(-1, 1),
      expand = expansion(0)
    )

  tabl_isotype_reduction <-
    isotype_data %>%
    group_by(sample_alias) %>%
    summarise(
      median_isotype_reduction = median(isotype_reduction) * 100,
      median_pre_denoise_isotype = median(pre_denoise_isotype_fraction) * 100,
      median_post_denoise_isotype = median(isotype_fraction) * 100,
      .groups = "drop"
    ) %>%
    mutate(across(where(is.numeric), ~ round(., 3))) %>%
    select(
      `Sample ID` = sample_alias,
      `Median isotype reduction [%]` = median_isotype_reduction,
      `Median pre-denoise isotype [%]` = median_pre_denoise_isotype,
      `Median post-denoise isotype [%]` = median_post_denoise_isotype
    ) %>%
    style_table(caption = "Isotype reduction by denoising", interactive = FALSE)

  list(
    plots = list(
      removed_umis = p_removed_umis,
      by_method = p_by_method,
      isotype_reduction = p_isotype_reduction
    ),
    tables = list(
      removed_umis = tabl_removed_umis,
      by_method = tabl_by_method,
      isotype_reduction = tabl_isotype_reduction
    )
  )
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
  pg_data
) {
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

#' Create the component for molecule plot
#'
#' This function generates a violin plot showing the distribution of the number
#' of molecules per sample, with a horizontal line indicating the molecule count cutoff.
#'
#' @param pg_data A Seurat object containing the data to be plotted.
#' @param sample_palette A named vector of colors for the samples.
#'
#' @return A ggplot object representing the violin plot of molecule counts.
#'
#' @export
#'
component_molecule_plot <- function(
  pg_data,
  sample_palette
) {
  p <- pg_data[[]] %>%
    plot_violin(
      x = "sample_alias",
      y = "n_umi",
      y_label = "Number of molecules",
      round = 0,
      use_log10 = TRUE,
      palette = sample_palette,
      alpha = 1,
      title = "Protein molecule counts",
    )
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
    if ("qc_files" %in% names(sample_qc_metrics)) {
      qc_data <- if (is.null(sample_qc_metrics$pool_qc_files)) {
        sample_qc_metrics$qc_files
      } else {
        sample_qc_metrics$pool_qc_files
      }
    } else {
      qc_data <- sample_qc_metrics
    }

    plot_data <-
      qc_data %>%
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

    plot_data <- set_sample_levels(plot_data, sample_levels)

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
        y_label = "Reads per cell (After filtering)",
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


#' Create the component sequencing saturation curve
#'
#' This function creates a component that visualizes the sequencing saturation
#' curve for each sample.
#'
#' @param object A Seurat object containing the sample data.
#' @param data_files A data frame containing the sample aliases and corresponding file names.
#' @param seqsat_comps The number of components to consider for sequencing saturation.
#' @param sample_levels Optional vector of sample levels to order the samples in the plots.
#'
#' @return A list containing plots showing the sequencing saturation curve for
#' each sample.
#'
#' @export
#'
component_sequencing_saturation_curve <-
  function(object, data_files, seqsat_comps = 100L, sample_levels = NULL) {
    object_meta <-
      FetchData(object, c("sample_alias", "reads_in_component")) %>%
      group_by(sample_alias) %>%
      summarise(
        n_comps = n(),
        reads_per_component = mean(reads_in_component)
      )

    sample_frac <- rev(seq(0.1, 1, 0.1))

    seqsat_curve_data <-
      data_files$sample_alias %>%
      set_names() %>%
      lapply(function(sampl) {
        n_comps <-
          object_meta$n_comps[which(object_meta$sample_alias == sampl)]
        n_comps <- pmin(n_comps, seqsat_comps)
        sampled_comps <- sample(
          object[[]] %>%
            filter(sample_alias == sampl) %>%
            rownames() %>%
            sample(n_comps)
        ) %>%
          # Strip the sample alias prefix from the component names
          # to match the component names in the pxl file database
          stringr::str_remove(paste0(sampl, "_"))

        pxl_file <- data_files %>%
          filter(sample_alias == sampl) %>%
          pull(filename)
        db <- PixelDB$new(pxl_file)
        on.exit(db$close())
        approximate_saturation_curve(
          db,
          fracs = sample_frac,
          components = sampled_comps,
          node_reads_multiplier = 1,
          verbose = FALSE
        )
      }) %>%
      bind_rows(.id = "sample_alias") %>%
      left_join(
        object_meta,
        by = join_by(sample_alias)
      ) %>%
      mutate(
        reads_per_component = reads_per_component * p
      ) %>%
      select(
        sample_alias,
        reads_per_component,
        graph_node_saturation = node_saturation,
        graph_edge_saturation = edge_saturation
      ) %>%
      pivot_longer(
        cols = c("graph_node_saturation", "graph_edge_saturation"),
        names_to = "type", values_to = "saturation"
      ) %>%
      mutate(type = str_remove(type, "graph_") %>%
        str_remove("_saturation") %>%
        str_to_sentence())

    seqsat_curve_data <-
      set_sample_levels(
        seqsat_curve_data,
        sample_levels
      ) %>%
      mutate(saturation = ifelse(
        saturation < 0, NA_real_, saturation
      )) %>%
      mutate(saturation = saturation * 100)

    seqsat_curve_data_mean <-
      seqsat_curve_data %>%
      group_by(sample_alias, type, reads_per_component) %>%
      summarise(
        saturation_mean = mean(saturation, na.rm = TRUE),
        saturation_se = sd(saturation, na.rm = TRUE) / sqrt(n()),
        .groups = "drop"
      ) %>%
      ungroup()

    plots <-
      seqsat_curve_data_mean %>%
      pull(sample_alias) %>%
      levels() %>%
      set_names() %>%
      lapply(function(x) {
        df_bg <- seqsat_curve_data_mean %>%
          filter(sample_alias != x) %>%
          stats::na.omit()
        df_fg <- seqsat_curve_data_mean %>%
          filter(sample_alias == x) %>%
          stats::na.omit()

        seqsat_curve_data_mean %>%
          ggplot() +
          # Background curves for other samples
          geom_line(
            data = df_bg,
            aes(reads_per_component, saturation_mean, group = sample_alias),
            color = "gray80"
          ) +
          # Foreground curve for the current sample
          geom_errorbar(
            data = df_fg,
            aes(
              x = reads_per_component,
              ymin = saturation_mean - saturation_se,
              ymax = saturation_mean + saturation_se
            ),
            color = "black",
            width = 0
          ) +
          geom_line(
            data = df_fg,
            aes(reads_per_component, saturation_mean),
            color = "black"
          ) +
          geom_point(
            data = df_fg,
            aes(reads_per_component, saturation_mean),
            color = "black",
            size = 0.5
          ) +
          facet_wrap(~type) +
          scale_y_continuous(limits = c(0, 100)) +
          theme_bw() +
          theme(
            axis.text.x = element_text(angle = 60, hjust = 1),
            panel.grid = element_blank(),
            legend.position = "none"
          ) +
          labs(
            x = "Reads per component",
            y = "Sequencing saturation [%]",
            title = "Sequencing saturation curve"
          )
      })

    return(plots)
  }

#' Create the component cell recovery
#'
#' This function creates a component that visualizes the cell recovery
#' after filtering components based on size.
#'
#' @param object A Seurat object containing the sample data.
#' @param sample_qc_metrics A list of sample QC metrics.
#' @param sample_levels Optional vector of sample levels to order the samples in the plots.
#'
#' @return A list containing a plot and a table summarizing the cell recovery.
#'
#' @export
#'
component_cell_recovery <-
  function(object, sample_qc_metrics, sample_levels = NULL) {
    if (is.null(sample_qc_metrics$pool_qc_files)) {
      qc_data <- sample_qc_metrics$qc_files
      hashed_experiment <- FALSE
    } else {
      qc_data <- sample_qc_metrics$pool_qc_files
      hashed_experiment <- TRUE
    }

    # Number of components before and after filtering, and size thresholds

    plot_data1 <-
      extract_sample_qc_metrics(
        qc_data,
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

    plot_data1 <- set_sample_levels(plot_data1, sample_levels)

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
      mutate(
        type = str_extract(type, "pre|post") %>%
          str_to_sentence() %>%
          paste("filtering") %>%
          factor(c("Pre filtering", "Post filtering")),
        label = ifelse(n > 1e5,
          compact_num(n),
          as.character(n)
        )
      ) %>%
      ggplot(aes(sample_alias, n)) +
      geom_col(fill = "#DAD6D7") +
      geom_text(aes(label = label),
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


    # Molecule Rank Plots

    plot_data_thresholds <-
      extract_sample_qc_metrics(
        qc_data,
        c(
          min_size_theshold =
            "component_size_min_filtering_threshold",
          max_size_theshold =
            "component_size_max_filtering_threshold"
        ),
        "graph"
      )

    plot_data <-
      qc_data %>%
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

    plot_data <- set_sample_levels(plot_data, sample_levels)

    .plot_MRP <-
      function(plot_data, plot_data_thresholds, selected_sample) {
        thresh <- plot_data_thresholds %>%
          filter(sample_alias == selected_sample)

        ggplot(plot_data, aes(rank, nodes, color = sample_alias == selected_sample)) +
          geom_point(size = 0.1, show.legend = FALSE) +
          geom_hline(
            data = thresh,
            aes(yintercept = min_size_theshold),
            linetype = "dashed"
          ) +
          geom_text(
            data = thresh,
            aes(
              x = 1,
              y = min_size_theshold,
              label = min_size_theshold
            ),
            vjust = -0.5,
            hjust = 0
          ) +
          geom_hline(
            data = thresh,
            aes(yintercept = max_size_theshold),
            linetype = "dashed"
          ) +
          geom_text(
            data = thresh,
            aes(
              x = 1,
              y = max_size_theshold,
              label = max_size_theshold
            ),
            vjust = -0.5,
            hjust = 0
          ) +
          scale_x_log10() +
          scale_y_log10() +
          scale_color_manual(values = c("TRUE" = "black", "FALSE" = "gray80")) +
          theme_bw() +
          theme(legend.position = "none") +
          labs(
            x = "Component rank (by number of molecules)",
            y = "Number of molecules",
            title = "Molecule rank plot"
          )
      }

    MRP_plots <-
      plot_data %>%
      pull(sample_alias) %>%
      levels() %>%
      set_names() %>%
      lapply(function(x) {
        plot_data %>%
          arrange(sample_alias == x) %>%
          .plot_MRP(plot_data_thresholds, x)
      })

    if (hashed_experiment) {
      plot_data_sample <-
        object[[]] %>%
        select(sample_alias, nodes = n_umi) %>%
        arrange(sample_alias, -nodes) %>%
        group_by(sample_alias) %>%
        mutate(rank = row_number()) %>%
        ungroup()

      plot_data_sample <- set_sample_levels(plot_data_sample, sample_levels)

      plots_sample <-
        plot_data_sample %>%
        pull(sample_alias) %>%
        levels() %>%
        set_names() %>%
        lapply(function(x) {
          plot_data_sample %>%
            arrange(sample_alias == x) %>%
            .plot_MRP(plot_data_thresholds, x)
        })

      MRP_plots <-
        list(
          "Pools" = MRP_plots,
          "Samples" = plots_sample
        )
    }

    # Table

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
      plots = list(p1, MRP_plots),
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

    plot_data <- set_sample_levels(plot_data, sample_levels)

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

    plot_data <- set_sample_levels(plot_data, sample_levels)

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
      rename(sample_alias = 1) %>%
      mutate(label = paste0(round(percent, 2), "%"))

    plot_data <- set_sample_levels(plot_data, sample_levels)

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
            round(100 * removed_total / total_edges_in, 2), "%"
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
      rename(sample_alias = 1) %>%
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
            select(sample_component, molecules) %>%
            distinct() %>%
            arrange(-molecules, sample_component) %>%
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

#' Create the component for abundance per marker
#'
#' This function creates plots visualizing the abundance of each marker
#' across samples and, optionally, across cell types.
#'
#' @param object A processed data object containing marker abundance data.
#' @param params A list of parameters, including control markers.
#' @param sample_palette A color palette for the samples.
#' @param per_celltype A boolean indicating whether to add a plot per cell type
#'   (default is TRUE).
#' @param sample_levels Optional vector of sample levels to order the samples in the plots.
#' @param test_mode A boolean indicating whether to run in test mode (default is FALSE).
#'
#' @return A list containing plots for each marker.
#'
#' @export
#'
component_abundance_per_marker <- function(
  object,
  params,
  sample_palette,
  per_celltype = TRUE,
  sample_levels = NULL,
  test_mode = FALSE
) {
  .abundance_plot_data <- function(object, params, test_mode = FALSE) {
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
            "condition",
            "seurat_clusters",
            "celltype"
          )
        ) %>%
          as_tibble(rownames = "cell_id"),
        by = "cell_id"
      ) %>%
      mutate(
        marker = factor(marker, order_cd_markers(
          unique(marker),
          params$control_markers
        ))
      ) %>%
      group_by(marker)
  }

  plot_data <- .abundance_plot_data(object, params, test_mode)

  plots <-
    plot_data %>%
    group_split() %>%
    set_names(group_keys(plot_data)$marker) %>%
    lapply(function(g_data) {
      p1 <-
        g_data %>%
        set_sample_levels(sample_levels) %>%
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
          use_grid = TRUE
        ) +
        guides(fill = "none")

      if (!per_celltype) {
        return(p1)
      }

      p2 <-
        g_data %>%
        set_sample_levels(sample_levels) %>%
        filter(celltype %in% displayed_cell_types) %>%
        complete(
          sample_alias = levels(sample_alias),
          celltype = displayed_cell_types
        ) %>%
        mutate(celltype = factor(celltype, displayed_cell_types)) %>%
        plot_violin(
          x = "sample_alias",
          y = "normcount",
          fill = "condition",
          y_label = "Normalized counts",
          summarize = FALSE,
          palette = sample_palette,
          facet_var = "celltype",
          use_jitter = TRUE,
          use_grid = TRUE
        ) +
        guides(fill = "none")

      p1 / p2
    })

  return(plots)
}


#' Create the component for proximity Z scores
#'
#' This function creates plots visualizing the proximity Z scores for selected
#' contrasts.
#'
#' @param object A Seurat object containing the sample data.
#' @param proximity_scores A data frame containing proximity scores for different markers.
#' @param sample_palette A color palette for the samples.
#' @param proximity_score One of "join_count_z" or "log2_ratio".
#' @param selected_contrasts A boolean indicating whether to filter for selected contrasts (default is TRUE).
#' @param sample_levels Optional vector of sample levels to order the samples in the plots.
#' @param test_mode A boolean indicating whether to run in test mode (default is FALSE).
#'
#' @return A list containing plots.
#'
#' @export
#'
component_proximity_selected <-
  function(
    object,
    proximity_scores,
    sample_palette,
    proximity_score = "log2_ratio",
    selected_contrasts = TRUE,
    sample_levels = NULL,
    test_mode = FALSE
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
      enframe("marker_1", "marker_2")

    processed_data <-
      proximity_scores %>%
      filter(celltype %in% displayed_cell_types)

    if (selected_contrasts) {
      processed_data <-
        processed_data %>%
        inner_join(plot_contrasts)
    }

    processed_data <-
      processed_data %>%
      set_sample_levels(sample_levels = sample_levels) %>%
      unite("contrast", marker_1, marker_2, sep = "/", remove = FALSE) %>%
      group_by(contrast)

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
            celltype = displayed_cell_types,
            fill = setNames(list(NA), proximity_score)
          ) %>%
          ggplot(aes(x = sample_alias, y = !!sym(proximity_score))) +
          geom_hline(yintercept = 0) +
          geom_violin(aes(fill = condition),
            draw_quantiles = 0.5,
            drop = FALSE,
            scale = "width"
          ) +
          geom_jitter(
            size = 0.1,
            position = position_jitter(height = 0),
            alpha = 0.3
          ) +
          geom_text(
            data = . %>%
              group_by(sample_alias, celltype) %>%
              summarise(!!sym(proximity_score) := median(!!sym(proximity_score)), .groups = "drop"),
            aes(label = round(!!sym(proximity_score), 2)),
            size = 3,
            vjust = 0
          ) +
          facet_grid(celltype ~ .) +
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
#' @param proximity_scores A data frame containing proximity scores for different markers.
#' @param sample_palette A color palette for the samples.
#' @param per_celltype A boolean indicating whether to add a plot per cell type (default is TRUE).
#' @param sample_levels Optional vector of sample levels to order the samples in the plots.
#' @param test_mode A boolean indicating whether to run in test mode (default is FALSE).
#'
#' @return A list containing plots for each marker.
#'
#' @export
#'
component_proximity_per_marker <- function(
  proximity_scores,
  sample_palette,
  per_celltype = TRUE,
  sample_levels = NULL,
  test_mode = FALSE
) {
  plot_data <-
    proximity_scores %>%
    filter(marker_1 == marker_2) %>%
    {
      if (test_mode) {
        filter(., marker_1 %in% head(levels(marker_1), 10))
      } else {
        .
      }
    } %>%
    mutate(
      celltype = factor(celltype, displayed_cell_types)
    ) %>%
    group_by(marker_1)

  plots <-
    plot_data %>%
    group_split() %>%
    set_names(group_keys(plot_data)$marker_1) %>%
    lapply(function(g_data) {
      p1 <-
        g_data %>%
        complete(
          sample_alias = levels(g_data$sample_alias)
        ) %>%
        set_sample_levels(sample_levels = sample_levels) %>%
        plot_violin(
          x = "sample_alias",
          y = "log2_ratio",
          fill = "condition",
          title = unique(g_data$marker_1) %>% as.character(),
          y_label = "Log2 ratio Proximity Score",
          summarize = FALSE,
          palette = sample_palette,
          use_jitter = TRUE,
          use_grid = TRUE,
          hline = 0
        ) +
        guides(fill = "none")

      if (!per_celltype) {
        return(p1)
      }

      p2 <-
        g_data %>%
        filter(celltype %in% displayed_cell_types) %>%
        complete(
          sample_alias = levels(g_data$sample_alias),
          celltype = displayed_cell_types
        ) %>%
        set_sample_levels(sample_levels = sample_levels) %>%
        plot_violin(
          x = "sample_alias",
          y = "log2_ratio",
          fill = "condition",
          y_label = "Log2 ratio Proximity Score",
          summarize = FALSE,
          palette = sample_palette,
          facet_var = "celltype",
          use_jitter = TRUE,
          use_grid = TRUE,
          hline = 0
        ) +
        guides(fill = "none")

      p1 / p2
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
        "celltype",
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
      sample_alias, condition, marker, celltype,
      mean_join_count_z, mean_log2_ratio
    )

  plot_join_count_z <-
    plot_data %>%
    ggplot(aes(marker, celltype, fill = mean_join_count_z)) +
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
    ggplot(aes(marker, celltype, fill = mean_log2_ratio)) +
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
#' @param proximity_scores A data frame containing proximity scores for different markers.
#' @param pg_data_processed A processed data object.
#' @param heatmap_gradient A color palette for the heatmaps.
#' @param min_pct_detected Minimum percentage of cells in which a marker must be detected to be included
#' (default is 0.25).
#' @param n_markers The number of markers to plot (default is 40).
#' @param plot_markers A vector of markers to plot (default is NULL).
#' @param test_mode A boolean indicating whether to run in test mode (default is FALSE).
#'
#' @return A list containing plots.
#'
#' @export
#'
component_proximity_heatmap_sample <- function(
  proximity_scores,
  pg_data_processed,
  heatmap_gradient,
  n_markers = 40,
  min_pct_detected = 0.25,
  plot_markers = NULL,
  test_mode = FALSE
) {
  pixelatorR:::assert_class(proximity_scores, "tbl_df")
  pixelatorR:::assert_class(pg_data_processed, "Seurat")
  pixelatorR:::assert_class(heatmap_gradient, "character")
  pixelatorR:::assert_class(n_markers, "numeric")
  pixelatorR:::assert_within_limits(n_markers, c(1, 100))
  pixelatorR:::assert_class(plot_markers, "character", allow_null = TRUE)
  pixelatorR:::assert_class(test_mode, "logical")

  processed_data <-
    summarize_colocalization_scores_per_sample(
      proximity_scores,
      plot_markers = plot_markers,
      test_mode = test_mode
    ) %>%
    bind_rows(select(.,
      marker_2 = marker_1,
      marker_1 = marker_2,
      everything()
    )) %>%
    distinct()

  if (is.null(plot_markers)) {
    plot_markers <-
      find_top_abundance_markers(pg_data_processed,
        n_markers = n_markers,
        summary_method = "mean"
      ) %>%
      pull(marker)
  }

  # Filter and symmetrise data
  processed_data <-
    processed_data %>%
    filter(marker_1 %in% plot_markers & marker_2 %in% plot_markers) %>%
    mutate(pct_detected = 100 * pct_detected) %>%
    group_by(sample_alias)


  plot_order <-
    processed_data %>%
    group_by(marker_1, marker_2) %>%
    summarise(
      mean_log2_ratio = mean(mean_log2_ratio, na.rm = TRUE)
    ) %>%
    ungroup() %>%
    pivot_wider(
      names_from = marker_2,
      values_from = mean_log2_ratio,
      values_fill = 0
    ) %>%
    column_to_rownames("marker_1") %>%
    as.matrix() %>%
    dist() %>%
    hclust(method = "ward.D2") %>%
    with(labels[order])

  plots <-
    processed_data %>%
    group_split() %>%
    set_names(group_keys(processed_data)$sample_alias) %>%
    lapply(function(sample_data) {
      samp_id <- unique(sample_data$sample_alias)
      cond <- unique(sample_data$condition)

      sample_data %>%
        mutate(
          marker_1 = factor(marker_1, levels = plot_order),
          marker_2 = factor(marker_2, levels = plot_order)
        ) %>%
        arrange(marker_1, marker_2) %>%
        ColocalizationHeatmap(
          marker1_col = "marker_1",
          marker2_col = "marker_2",
          value_col = "mean_log2_ratio",
          type = "dots",
          size_col = "pct_detected",
          cluster_rows = FALSE,
          cluster_cols = FALSE,
          symmetrise = FALSE,
          cellwidth = 8,
          cellheight = 8,
          fontsize = 8,
          colors = heatmap_gradient,
          border_color = "white",
          legend_range = c(-1, 1) * 0.75
        ) +
        scale_size_continuous(range = c(0, 5), limits = c(0, 100)) +
        labs(
          title = paste(samp_id),
          subtitle = cond,
          fill = "Mean\nProximity\nlogratio",
          size = "% cells\ndetected"
        ) +
        theme(axis.title = element_blank())
    })

  return(plots)
}


#' Create the component for proximity heatmaps per celltype and condition
#'
#' This function creates plots visualizing the proximity scores as heatmaps per celltype and condition.
#'
#' @param proximity_scores A data frame containing proximity scores for different markers.
#' @param pg_data_processed A processed data object.
#' @param heatmap_gradient A color palette for the heatmaps.
#' @param n_markers The number of markers to plot (default is 40).
#' @param min_pct_detected Minimum percentage of cells in which a marker must be detected to be included
#' (default is 0.25).
#' @param plot_markers A vector of markers to plot (default is NULL).
#' @param test_mode A boolean indicating whether to run in test mode (default is FALSE).
#'
#' @return A list containing plots.
#'
#' @export
#'
component_proximity_heatmap_celltype <- function(
  proximity_scores,
  pg_data_processed,
  heatmap_gradient,
  n_markers = 40,
  min_pct_detected = 0.25,
  plot_markers = NULL,
  test_mode = FALSE
) {
  pixelatorR:::assert_class(proximity_scores, "tbl_df")
  pixelatorR:::assert_class(pg_data_processed, "Seurat")
  pixelatorR:::assert_class(heatmap_gradient, "character")
  pixelatorR:::assert_class(n_markers, "numeric")
  pixelatorR:::assert_within_limits(n_markers, c(1, 100))
  pixelatorR:::assert_class(plot_markers, "character", allow_null = TRUE)
  pixelatorR:::assert_class(test_mode, "logical")

  processed_data <-
    proximity_scores %>%
    filter(celltype %in% displayed_cell_types) %>%
    summarize_colocalization_scores_per_celltype(
      plot_markers = plot_markers,
      test_mode = test_mode
    ) %>%
    bind_rows(select(.,
      marker_2 = marker_1,
      marker_1 = marker_2,
      everything()
    )) %>%
    distinct()

  if (is.null(plot_markers)) {
    top_markers <-
      find_top_abundance_markers(pg_data_processed,
        n_markers = n_markers,
        summary_method = "mean",
        group_col = "celltype"
      )

    processed_data <-
      processed_data %>%
      inner_join(top_markers, by = c("celltype", "marker_1" = "marker")) %>%
      inner_join(top_markers, by = c("celltype", "marker_2" = "marker"))
  } else {
    processed_data <-
      processed_data %>%
      filter(marker_1 %in% plot_markers & marker_2 %in% plot_markers)
  }

  processed_data <-
    processed_data %>%
    mutate(
      celltype = factor(celltype, displayed_cell_types),
      pct_detected = 100 * pct_detected
    ) %>%
    group_by(celltype) %>%
    arrange(celltype)

  plots <-
    processed_data %>%
    group_split() %>%
    set_names(group_keys(processed_data)$celltype) %>%
    lapply(function(cell_type_data) {
      cell_type_data <-
        cell_type_data %>%
        group_by(sample_alias)

      plot_order <-
        cell_type_data %>%
        group_by(marker_1, marker_2) %>%
        summarise(
          mean_log2_ratio = mean(mean_log2_ratio, na.rm = TRUE),
          .groups = "drop"
        ) %>%
        ungroup() %>%
        pivot_wider(
          names_from = marker_2,
          values_from = mean_log2_ratio,
          values_fill = 0
        ) %>%
        column_to_rownames("marker_1") %>%
        as.matrix() %>%
        dist() %>%
        hclust(method = "ward.D2") %>%
        with(labels[order])

      cell_type_data %>%
        group_split() %>%
        set_names(group_keys(cell_type_data)$sample_alias) %>%
        lapply(function(sample_data) {
          cellpop <- unique(sample_data$celltype)
          samp_id <- unique(sample_data$sample_alias)
          cond <- unique(sample_data$condition)

          sample_data %>%
            mutate(
              marker_1 = factor(marker_1, levels = plot_order),
              marker_2 = factor(marker_2, levels = plot_order)
            ) %>%
            arrange(marker_1, marker_2) %>%
            ColocalizationHeatmap(
              marker1_col = "marker_1",
              marker2_col = "marker_2",
              value_col = "mean_log2_ratio",
              type = "dots",
              size_col = "pct_detected",
              cluster_rows = FALSE,
              cluster_cols = FALSE,
              symmetrise = FALSE,
              cellwidth = 8,
              cellheight = 8,
              fontsize = 8,
              colors = heatmap_gradient,
              border_color = "white",
              legend_range = c(-1, 1) * 0.8
            ) +
            scale_size_continuous(range = c(0, 5), limits = c(0, 100)) +
            labs(
              title = paste(samp_id, "-", cellpop),
              subtitle = cond,
              fill = "Mean\nProximity\nlogratio",
              size = "% cells\ndetected"
            ) +
            theme(axis.title = element_blank())
        })
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
#' @param plot_height Numeric value specifying the height of the full plot. Used to space the legend.
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
    sample_plots = FALSE,
    plot_height = 5
  ) {
    plot_reductions %>%
      set_names(toupper(plot_reductions)) %>%
      lapply(function(red) {
        p1 <-
          plot_embedding(object, red,
            metavars = "seurat_clusters",
            pal = cluster_palette,
            plot_title = "Clusters"
          ) +
          theme(plot.title = element_text(hjust = 0.5))
        p2 <-
          plot_embedding(object, red,
            metavars = "condition",
            pal = sample_palette,
            plot_title = "Conditions",
            label = FALSE,
            legend_position = "bottom",
            yaxis_title = FALSE,
            extract_legend = TRUE,
            plot_height = plot_height / 4
          )

        p2$plot <-
          p2$plot +
          theme(plot.title = element_text(hjust = 0.5))

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
              Combined =
                p1 + p2$plot +
                  plot_void() + p2$legend +
                  plot_layout(nrow = 2, heights = c(3, 1)),
              Samplewise = samplewise_plots
            )
        } else {
          plot_list <-
            p1 + p2$plot +
            plot_void() + p2$legend +
            plot_layout(nrow = 2, heights = c(3, 1))
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
        metavars = "celltype",
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
        vars = c("seurat_clusters", "celltype")
      ) %>%
        as_tibble(rownames = "cell_id")) %>%
      group_by(seurat_clusters, cell_annotation = celltype, marker) %>%
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
          "celltype",
          "condition",
          "sample_alias"
        )
      ) %>%
      as_tibble(rownames = "cell_id") %>%
      group_by(sample_alias, condition, celltype) %>%
      count() %>%
      group_by(sample_alias, condition) %>%
      mutate(frac = n / sum(n)) %>%
      ungroup() %>%
      complete(
        nesting(sample_alias, condition),
        celltype,
        fill = list(n = 0, frac = 0)
      ) %>%
      mutate(celltype = factor(celltype, displayed_cell_types)) %>%
      group_by(sample_alias)

    p4 <-
      plot_data %>%
      ggplot(aes(sample_alias, frac, fill = celltype)) +
      geom_col(position = "stack") +
      geom_text(aes(label = scales::percent(frac, accuracy = 0.1)),
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
          ggplot(aes(celltype, n, fill = celltype)) +
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
      arrange(celltype) %>%
      group_by(celltype)

    barplots2 <-
      plot_data %>%
      group_split() %>%
      set_names(group_keys(plot_data)$celltype) %>%
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
              unique(g_data$celltype),
              "cells per sample"
            )
          )
      })

    tabl1 <-
      plot_data %>%
      ungroup() %>%
      mutate(frac = scales::percent(frac, accuracy = 0.1)) %>%
      select(
        `Sample ID` = sample_alias,
        `Cell annotation` = celltype,
        frac
      ) %>%
      pivot_wider(names_from = `Cell annotation`, values_from = frac, values_fill = "0.0%") %>%
      style_table(caption = "Cell type composition [%]", interactive = FALSE)

    tabl2 <-
      plot_data %>%
      ungroup() %>%
      select(
        `Sample ID` = sample_alias,
        `Cell annotation` = celltype,
        n
      ) %>%
      pivot_wider(names_from = `Cell annotation`, values_from = n, values_fill = 0) %>%
      style_table(caption = "Cell type composition", interactive = FALSE)

    return(list(
      dimred_plots = list(
        "Level 1 annotation" = p1,
        "Level 2 annotation" = p2
      ),
      heatmap = p3,
      celltype_composition = p4,
      celltype_composition_barplots1 = barplots1,
      celltype_composition_barplots2 = barplots2,
      celltype_composition_table = tabl1,
      celltype_numbers_table = tabl2
    ))
  }


#' Create the component sequencing saturation
#'
#' This function creates a component that visualizes the sequencing saturation
#' for each sample.
#'
#' @param qc_metrics_tables A list of QC metrics tables.
#' @param sample_levels Optional vector of sample levels to order the samples in the plots.
#'
#' @return A list containing plots and a table summarizing the sequencing saturation.
#'
#' @export
#'
component_sequencing_saturation <-
  function(qc_metrics_tables, sample_levels = NULL) {
    ss <- qc_metrics_tables$seq_saturation %>%
      rename(sample_alias = 1)
    ss <- set_sample_levels(ss, sample_levels)

    p1 <-
      ss %>%
      ggplot(aes(sample_alias, fraction_valid_reads)) +
      geom_col(fill = "#DAD6D7") +
      geom_text(aes(label = paste0(round(fraction_valid_reads, 3), " %")),
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
        x = NULL, y = "Valid reads [%]",
        title = "Fraction of valid reads",
        subtitle = "Valid reads from total reads",
        caption = expression(~ italic(frac("# Valid reads", "# Total reads")))
      )

    p2 <-
      ss %>%
      ggplot(aes(sample_alias, fraction_graph_reads)) +
      geom_col(fill = "#DAD6D7") +
      geom_text(aes(label = paste0(round(fraction_graph_reads, 3), " %")),
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
        x = NULL, y = "Graph reads [%]",
        title = "Fraction of graph reads",
        subtitle = "Reads in valid components from total reads",
        caption = expression(~ italic(frac("# Graph reads", "# Total reads")))
      )

    p3 <-
      ss %>%
      ggplot(aes(sample_alias, valid_reads_saturation)) +
      geom_col(fill = "#DAD6D7") +
      geom_text(aes(label = paste0(round(valid_reads_saturation, 3), " %")),
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
        x = NULL, y = "Sequencing saturation [%]",
        title = "Valid read saturation",
        subtitle = "Redundancy in valid reads",
        caption = expression(~ italic("Saturation" == 1 -
          frac(
            "# Deduped valid reads",
            "# Valid reads"
          )))
      )

    p4 <-
      ss %>%
      ggplot(aes(sample_alias, graph_edge_saturation)) +
      geom_col(fill = "#DAD6D7") +
      geom_text(aes(label = paste0(round(graph_edge_saturation, 3), " %")),
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
        x = NULL, y = "Sequencing saturation [%]",
        title = "Graph edge saturation",
        subtitle = "Redundancy in graph reads",
        caption = expression(~ italic("Saturation" == 1 -
          frac(
            "# Graph edges",
            "# Graph reads"
          )))
      )

    p5 <-
      ss %>%
      ggplot(aes(sample_alias, graph_node_saturation)) +
      geom_col(fill = "#DAD6D7") +
      geom_text(aes(label = paste0(round(graph_node_saturation, 3), " %")),
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
        x = NULL, y = "Sequencing saturation [%]",
        title = "Graph node saturation",
        subtitle = "Redundancy in graph proteins",
        caption = expression(~ italic("Saturation" == 1 -
          frac(
            "# Graph proteins",
            "# Graph reads"
          )))
      )

    tabl <-
      ss %>%
      arrange(sample_alias) %>%
      mutate(across(c(
        valid_reads_saturation, graph_node_saturation,
        graph_edge_saturation, fraction_valid_reads,
        fraction_graph_reads
      ), ~ round(., 3))) %>%
      mutate(across(
        c(
          graph_proteins, total_reads, valid_reads,
          deduped_valid_reads, graph_edges
        ),
        . %>%
          `/`(1e6) %>%
          round(3)
      )) %>%
      select(
        `Sample ID` = sample_alias,
        `Protein molecules [M]` = graph_proteins,
        `Total reads [M]` = total_reads,
        `Valid reads [M]` = valid_reads,
        `Deduped valid reads [M]` = deduped_valid_reads,
        `Total edges [M]` = graph_edges,
        `Fraction valid reads [%]` = fraction_valid_reads,
        `Fraction graph reads [%]` = fraction_graph_reads,
        `Valid reads saturation [%]` = valid_reads_saturation,
        `Graph edge saturation [%]` = graph_edge_saturation,
        `Graph node saturation [%]` = graph_node_saturation
      ) %>%
      style_table(caption = "Sequencing saturation", interactive = FALSE)

    return(list(
      plots = list(
        p1, p2, p3, p4, p5
      ),
      table = tabl
    ))
  }


#' Create the component hashing
#'
#' @param qc_metrics_tables QC metrics from [get_qc_metrics()] including `sample_hash_stats`.
#' @param colors A gradient color palette to use.
#'
#' @return List with `plot`, `table`, `heatmap_plots_hash_purity`, and `heatmap_plots_hash_fraction`.
#'
#' @export
#'
component_hashing <-
  function(qc_metrics_tables, colors) {
    component_stats <- qc_metrics_tables$sample_hash_stats$component_stats
    pixelatorR:::assert_class(component_stats, "tbl_df")
    n_unique_hashes <- component_stats$id %>%
      unique() %>%
      length()

    sample_stats <- qc_metrics_tables$sample_hash_stats$sample_stats
    pixelatorR:::assert_class(sample_stats, "tbl_df")

    component_stats_heatmap_purity <- qc_metrics_tables$sample_hash_stats$component_stats_heatmap_purity
    pixelatorR:::assert_class(component_stats_heatmap_purity, "tbl_df")

    component_stats_heatmap_fraction <- qc_metrics_tables$sample_hash_stats$component_stats_heatmap_fraction
    pixelatorR:::assert_class(component_stats_heatmap_fraction, "tbl_df")

    id_hash_order <- attributes(component_stats_heatmap_purity)$id_hash_order

    # Sample determined component purity
    component_hash_purity <- component_stats %>%
      group_by(component, sample_alias, version) %>%
      summarize(count = sum(count), .groups = "keep") %>%
      ungroup(version) %>%
      mutate(hash_fraction = count / sum(count)) %>%
      group_by(component, sample_alias)

    sample_top_id <- component_hash_purity %>%
      group_by(sample_alias) %>%
      dplyr::slice_max(hash_fraction, n = 1, with_ties = FALSE) %>%
      ungroup() %>%
      select(sample_alias, version) %>%
      distinct()

    plot_data <-
      component_hash_purity %>%
      dplyr::slice_max(hash_fraction, n = 1, with_ties = FALSE)

    # Use 90% as lower limit unless there are smaller values
    plot_limits <- c(floor(10 * min(c(min(plot_data$hash_fraction), 0.9))) / 10, 1)
    p1 <- suppressMessages(
      plot_violin(
        plot_data,
        x = "sample_alias",
        y = "hash_fraction",
        use_pct = TRUE,
        round = 2,
        expand = c(0, 0.1),
        title = "Hash purity",
        subtitle = "Fraction of UMIs from the most abundant hash",
        y_label = "Hash purity"
      ) +
        scale_y_continuous(limits = plot_limits, labels = scales::percent)
    )

    # Sample confidence plot
    plot_data <-
      qc_metrics_tables$sample_hash_stats$component_sample_confidence %>%
      arrange(pool) %>%
      group_by(pool)

    .generate_confidence_plots <- function(plot_data, metric, yes_no_palette) {
      # 1. Define metric-specific plot settings
      if (metric == "hash_enrichment_factor") {
        rect_ymin <- 0
        rect_ymax <- 1
        text_y <- 1
        y_scale <- scale_y_log10()
        y_axis_label <- "Hash enrichment factor"
        plot_title <- "Component hash enrichment factor"
      } else if (metric == "sample_confidence") {
        rect_ymin <- -Inf
        rect_ymax <- 0
        text_y <- 0
        y_scale <- scale_y_continuous(limits = c(0, 1), labels = scales::percent)
        y_axis_label <- "Hash purity"
        plot_title <- "Component sample confidence"
      }

      # 2. Generate the split lists and plots
      plot_data %>%
        group_split() %>%
        set_names(group_keys(plot_data)$pool) %>%
        lapply(function(g_data) {
          # Use tidy evaluation (!!sym(metric)) to arrange dynamically
          g_data <- g_data %>%
            arrange(desc(!!sym(metric))) %>%
            mutate(
              rank = row_number(),
              type = ifelse(sample_alias == "undetermined",
                "Undetermined", "Sample assigned"
              )
            )

          g_data_sum <- g_data %>%
            group_by(type) %>%
            count() %>%
            ungroup() %>%
            mutate(
              cumsum_n = cumsum(n),
              x = cumsum_n - n,
              hjust = c(0, 1)[seq_len(n())],
              text_pos = range(c(x, cumsum_n))[seq_len(n())]
            )

          # Use tidy evaluation (!!sym(metric)) for the y-axis aesthetic
          g_data %>%
            ggplot(aes(x = rank, y = !!sym(metric), color = type)) +
            geom_point(size = 0.5) +
            geom_rect(
              data = g_data_sum,
              aes(xmin = x, xmax = cumsum_n, ymin = rect_ymin, ymax = rect_ymax, fill = type),
              inherit.aes = FALSE,
              color = NA
            ) +
            geom_text(
              data = g_data_sum,
              aes(
                x = text_pos, y = text_y, label = paste0(type, "\n", n, " cells"),
                hjust = hjust
              ),
              vjust = -0.5,
              inherit.aes = FALSE,
              size = 3
            ) +
            theme_bw() +
            scale_color_manual(values = yes_no_palette, name = "") +
            scale_fill_manual(values = yes_no_palette, name = "") +
            y_scale + # Apply the dynamic scale defined earlier
            labs(x = "Component rank", y = y_axis_label, title = plot_title) +
            theme(
              legend.position = "none",
              panel.grid = element_blank()
            )
        })
    }

    metric <- intersect(c("sample_confidence", "hash_enrichment_factor"), names(plot_data))
    if (length(metric) == 0) {
      cli_abort("component_sample_confidence must contain either 'sample_confidence' or 'hash_enrichment_factor'.")
    }

    p2 <- .generate_confidence_plots(
      plot_data,
      metric[[1]],
      yes_no_palette = yes_no_palette
    )

    tabl <-
      sample_stats %>%
      select(
        `Sample ID` = sample_alias,
        `Hash purity B2M [%]` = mean_purity_B2M,
        `Hash purity CD298 [%]` = mean_purity_CD298,
        `Hash purity CD98 [%]` = mean_purity_CD98,
        `Fraction hash UMIs` = hash_pct
      ) %>%
      mutate(across(where(is.numeric), ~ round(.x * 100, 2))) %>%
      style_table(caption = "Hash purity", interactive = FALSE)

    component_stats_heatmap_purity <-
      component_stats_heatmap_purity %>%
      group_by(sample_alias)

    heatmap_plots_hash_purity <-
      component_stats_heatmap_purity %>%
      group_by(sample_alias) %>%
      group_split() %>%
      set_names(group_keys(component_stats_heatmap_purity)$sample_alias) %>%
      lapply(function(g_data) {
        if (n_distinct(g_data$sample_component) > 2000) {
          g_data <-
            g_data %>%
            head(2000)
        }

        sample_top_hash_id <- sample_top_id %>%
          filter(sample_alias == g_data$sample_alias[1]) %>%
          pull(version)

        heatmap_matrix <- g_data %>%
          select(sample_component, matches(paste0("-", sample_top_hash_id, "$"))) %>%
          column_to_rownames("sample_component") %>%
          as.matrix() %>%
          t() %>%
          (\(x) x * 100)
        rownames(heatmap_matrix) <- stringr::str_remove(rownames(heatmap_matrix), paste0("-", sample_top_hash_id, "$"))
        heatmap_matrix %>%
          ComplexHeatmap::pheatmap(
            color = colors,
            breaks = seq(0, 100, 1),
            cellheight = 10,
            show_colnames = FALSE,
            cluster_rows = FALSE,
            clustering_method = "ward.D2",
            heatmap_legend_param = list(title = "% Hash\npurity")
          ) %>%
          as.ggplot()
      })

    component_stats_heatmap_fraction <-
      component_stats_heatmap_fraction %>%
      group_by(sample_alias)

    heatmap_plots_hash_fraction <-
      component_stats_heatmap_fraction %>%
      group_split() %>%
      set_names(group_keys(component_stats_heatmap_fraction)$sample_alias) %>%
      lapply(function(g_data) {
        if (n_distinct(g_data$sample_component) > 2000) {
          g_data <-
            g_data %>%
            head(2000)
        }

        heatmap_matrix <- g_data %>%
          select(sample_component, where(is.numeric)) %>%
          column_to_rownames("sample_component") %>%
          as.matrix() %>%
          t() %>%
          (\(x) x * 100)

        heatmap_matrix <- heatmap_matrix[id_hash_order, , drop = FALSE]

        heatmap_matrix %>%
          ComplexHeatmap::pheatmap(
            color = colors,
            breaks = seq(0, 100, 1),
            cellheight = 10,
            show_colnames = FALSE,
            cluster_rows = FALSE,
            clustering_method = "ward.D2",
            heatmap_legend_param = list(title = "% fraction\nhash UMIs")
          ) %>%
          as.ggplot()
      })

    list(
      plot = p1,
      sample_confidence_plots = p2,
      table = tabl,
      heatmap_plots_hash_purity = heatmap_plots_hash_purity,
      heatmap_plots_hash_fraction = heatmap_plots_hash_fraction
    )
  }
