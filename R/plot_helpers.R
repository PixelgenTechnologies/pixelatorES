#' Header wrap plots
#'
#' Wraps a list of plots with headers for each plot.
#'
#' @param plots List of plots to be wrapped.
#' @param level Integer indicating the header level for the headers (default is 2).
#'
#' @return A formatted list of plots with headers.
#'
#' @export
title_plotlist <- function(plots, level = 2) {
  pixelatorR:::assert_class(plots, "list")
  for (plot in plots) pixelatorR:::assert_class(plot, c("ggplot", "datatables"))

  nams <- names(plots)

  if (is.null(nams)) nams <- rep(" ", length(plots))

  # Loop through each plot and create a tab for each
  for (tab in seq_len(length(plots))) {
    # Generate a header for each plot in the tabset
    cat(paste0(strrep("#", level), " ", nams[tab], "\n\n"))

    # Print the plot
    msg <- try(print(plots[[tab]]))
    if (inherits(msg, "try-error")) {
      print(
        ggplot() +
          geom_blank()
      )
    }
    cat("\n\n")
  }
}

#' Tabset plots
#'
#' Tabset a list of plots at a given header level.
#'
#' @param plots List of plots to be tabsetted.
#' @param level Integer indicating the header level for the tabset title (default is 2).
#'
#' @return A formatted tabset of plots.
#'
#' @export
tabset_plotlist <- function(plots, level = 2) {
  pixelatorR:::assert_class(plots, "list")
  for (plot in plots) pixelatorR:::assert_class(plot, c("ggplot", "datatables"))

  # Start the tabset
  cat("::: {.panel-tabset .nav-pills}\n")

  title_plotlist(plots, level)

  # Close the tabset
  cat(":::\n")
}


#' Tabset nested plots
#'
#' Tabset a list of plots, which may contain nested lists of plots, at a given header level.
#'
#' @param plots List of plots or list of lists of plots to be tabsetted.
#' @param level Integer indicating the header level for the tabset title (default is 2).
#'
#' @return A formatted tabset of plots.
#'
#' @export
tabset_nested_plotlist <- function(plots, level = 2) {
  pixelatorR:::assert_class(plots, "list")
  for (plot in plots) {
    pixelatorR:::assert_class(
      plot,
      c(
        "ggplot",
        "datatables",
        "list"
      )
    )
  }

  # Start the tabset
  cat("\n\n")
  cat("::: {.panel-tabset .nav-pills}\n")

  for (tab in seq_len(length(plots))) {
    if (inherits(plots[[tab]], "list")) {
      # If the plot is a list, call tabset_nested_plotlist recursively

      cat(paste0(strrep("#", level), " ", names(plots)[tab], "\n"))

      tabset_nested_plotlist(plots[[tab]], level + 1)
    } else if (inherits(plots[[tab]], "ggplot")) {
      # Otherwise, just print the plot
      list(plots[[tab]]) %>%
        set_names(names(plots)[tab]) %>%
        title_plotlist(level)
    }
  }
  # Close the tabset
  cat(":::\n")
}


#' Tabset a figure and a table
#'
#' Tabset a figure and a table at a given header level. After using this
#' function, you should close the tabset with `close_tabset()`.
#'
#' @param figure A ggplot object representing the figure to be tabsetted.
#' @param table A table to be tabsetted alongside the figure.
#' @param level Integer indicating the header level for the tabset title (default is 2).
#' @param mode Character indicating the mode of the tabset ("tabset" or "title").
#'
#' @return A formatted tabset containing the figure and table.
#'
#' @export
tabset_figure_table <- function(figure, table, level = 2,
                                  mode = "tabset") {
  pixelatorR:::assert_class(table, "datatables")
  mode <- match.arg(mode, c("tabset", "title"))

  set_func <-
    switch(mode,
      tabset = tabset_plotlist,
      title = title_plotlist
    )

  # Start the tabset
  cat("::: {.panel-tabset .nav-pills}\n")

  cat(paste0(strrep("#", level), " Figure\n\n"))

  if (inherits(figure, "list")) {
    for (plot in figure) pixelatorR:::assert_class(plot, "ggplot")

    set_func(figure, level + 1)
  } else {
    pixelatorR:::assert_class(figure, "ggplot")
    print(figure)
  }
  cat("\n\n")
  cat(paste0(strrep("#", level), " Table\n\n"))

  return(table)
}

#' Close Tabset
#'
#' Closes the tabset that was opened with `tabset_figure_table`.
#'
#' @return Nothing.
#'
#' @export
#'
close_tabset <-
  function() {
    cat(":::\n")
  }

#' Plot Embedding
#'
#' Plot a 2D embedding of a Seurat object, allowing for customization of the plot aesthetics.
#' @param object A Seurat object containing the embeddings to be plotted.
#' @param plot_reduction The name of the reduction to use for plotting (default is "pca").
#' @param dims A vector specifying the dimensions to plot (default is 1:2).
#' @param metavars A variable to color the points by (default is "seurat_clusters").
#' @param pal A color palette to use for the points (default is NULL).
#' @param label Logical indicating whether to label the points with their metadata (default is TRUE).
#' @param xaxis_title Logical indicating whether to show the x-axis title (default is TRUE).
#' @param yaxis_title Logical indicating whether to show the y-axis title (default is TRUE).
#' @param fix_aspect_ratio Logical indicating whether to fix the aspect ratio of the plot (default is TRUE).
#' @param plot_title The title of the plot (default is an empty string).
#' @param legend_position The position of the legend (default is "none").
#' @param legend_cols The number of columns in the legend (default is 2).
#'
#' @return A ggplot object representing the embedding.
#'
#' @export
plot_embedding <-
  function(object,
             plot_reduction = "pca",
             dims = 1:2,
             metavars = "seurat_clusters",
             pal = NULL,
             label = TRUE,
             plot_title = "",
             xaxis_title = TRUE,
             yaxis_title = TRUE,
             fix_aspect_ratio = TRUE,
             legend_position = "none",
             legend_cols = 2) {
    pixelatorR:::assert_class(object, "Seurat")
    pixelatorR:::assert_single_value(plot_reduction, type = "string")
    pixelatorR:::assert_vector(dims, "numeric", n = 2)
    pixelatorR:::assert_vector(metavars, "character", n = 1)
    pixelatorR:::assert_vector(pal, "character", n = 1, allow_null = TRUE)
    pixelatorR:::assert_single_value(label, type = "bool")
    pixelatorR:::assert_single_value(plot_title, type = "string")
    pixelatorR:::assert_single_value(legend_position, type = "string")
    pixelatorR:::assert_single_value(legend_cols, type = "numeric")

    embeddings <-
      Embeddings(object, reduction = plot_reduction) %>%
      as_tibble(rownames = "sample_component") %>%
      select(1, dims + 1)

    plot_metadata <-
      FetchData(object, vars = metavars) %>%
      as_tibble(rownames = "sample_component")

    plot_dimnames <-
      colnames(embeddings)[dims + 1]

    plot_data <-
      left_join(embeddings, plot_metadata,
        by = "sample_component"
      ) %>%
      rename(
        V1 = !!plot_dimnames[1], V2 = !!plot_dimnames[2],
        meta = !!metavars
      )


    if (xaxis_title) {
      plot_dimnames[1] <- beaut_label(plot_dimnames[1])
    } else {
      plot_dimnames[1] <- ""
    }

    if (yaxis_title) {
      plot_dimnames[2] <- beaut_label(plot_dimnames[2])
    } else {
      plot_dimnames[2] <- ""
    }

    p <-
      plot_data %>%
      ggplot(aes(x = V1, y = V2, color = meta)) +
      geom_point(size = 0.75) +
      theme_bw() +
      theme(
        panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = legend_position,
        legend.key.size = unit(0.5, "cm"),
        legend.spacing.x = unit(0.2, "cm"),
        legend.direction = "vertical"
      ) +
      coord_fixed() +
      labs(
        title = plot_title,
        x = plot_dimnames[1],
        y = plot_dimnames[2],
        color = ""
      ) +
      guides(color = guide_legend(
        ncol = legend_cols,
        override.aes = list(size = 3)
      ))

    if (!is.null(pal)) {
      p <- p + scale_color_manual(values = pal)
    }

    if (label) {
      p <-
        p +
        geom_text(
          data = . %>%
            group_by(meta) %>%
            summarise(V1 = median(V1), V2 = median(V2)),
          aes(label = meta),
          color = "black",
          alpha = 0.75
        )
    }

    if (fix_aspect_ratio) {
      V1_range <- range(plot_data$V1)
      V2_range <- range(plot_data$V2)

      range_diff <- diff(V2_range) - diff(V1_range)

      if (range_diff < 0) {
        V2_range <- V2_range + c(-1, 1) * abs(range_diff) / 2
      } else if (range_diff > 0) {
        V1_range <- V1_range + c(-1, 1) * abs(range_diff) / 2
      }

      p <-
        p +
        scale_x_continuous(limits = V1_range) +
        scale_y_continuous(limits = V2_range)
    }

    return(p)
  }

#' Plot Embeddings Samplewise
#'
#' Plot embeddings for each sample in a Seurat object, allowing for customization of the plot aesthetics.
#'
#' @param samples A vector of sample identifiers to plot.
#' @param object A Seurat object containing the embeddings to be plotted.
#' @param plot_reduction The name of the reduction to use for plotting (default is "pca").
#' @param dims A vector specifying the dimensions to plot (default is 1:2).
#' @param metavars A variable to color the points by (default is "seurat_clusters").
#' @param pal A color palette to use for the points (default is NULL).
#' @param label Logical indicating whether to label the points with their metadata (default is TRUE).
#' @param plot_title The title of the plot (default is an empty string).
#' @param legend_position The position of the legend (default is "none").
#' @param legend_cols The number of columns in the legend (default is 2).
#'
#' @return A list of ggplot objects, each representing the embedding for a sample.
#' @export
plot_embeddings_samplewise <-
  function(samples,
             object,
             plot_reduction = "pca",
             dims = 1:2,
             metavars = "seurat_clusters",
             pal = NULL,
             label = TRUE,
             plot_title = "",
             legend_position = "none",
             legend_cols = 2) {
    pixelatorR:::assert_vector(samples, "character", n = 1)
    pixelatorR:::assert_class(object, "Seurat")
    pixelatorR:::assert_single_value(plot_reduction, type = "string")
    pixelatorR:::assert_vector(dims, "numeric", n = 2)
    pixelatorR:::assert_vector(metavars, "character", n = 1)
    pixelatorR:::assert_vector(pal, "character", n = 1, allow_null = TRUE)
    pixelatorR:::assert_single_value(label, type = "bool")
    pixelatorR:::assert_single_value(plot_title, type = "string")
    pixelatorR:::assert_single_value(legend_position, type = "string")
    pixelatorR:::assert_single_value(legend_cols, type = "numeric")


    samples %>%
      set_names() %>%
      lapply(function(samp_id) {
        object %>%
          subset(subset = sample_alias == samp_id) %>%
          plot_embedding(
            plot_reduction = plot_reduction,
            dims = dims,
            metavars = metavars,
            pal = pal,
            label = label,
            plot_title = plot_title,
            legend_position = legend_position,
            legend_cols = legend_cols
          ) +
          scale_x_continuous(limits = range(Embeddings(object, plot_reduction)[, dims[1]])) +
          scale_y_continuous(limits = range(Embeddings(object, plot_reduction)[, dims[2]]))
      })
  }



#' Make a violin plot
#'
#' @param plot_data A data frame containing the data to plot.
#' @param x The variable to use for the x-axis.
#' @param y The variable to use for the y-axis.
#' @param fill Optional variable to fill color the violins.
#' @param palette Optional color palette to use for the fill variable.
#' @param title The title of the plot.
#' @param subtitle The subtitle of the plot.
#' @param x_label The label for the x-axis.
#' @param y_label The label for the y-axis.
#' @param fill_label The label for the fill variable.
#' @param summarize Logical indicating whether to show the median y value
#' above the 0.5 quantile.
#' @param round Number of decimal places to round the median y value.
#' @param use_pct Logical indicating whether to show the y-axis in percentage.
#' @param use_1k Logical indicating whether to divide median y value by 1,000.
#' @param expand Optional expansion factor for the y-axis.
#' @param use_log10 Logical indicating whether to use a log10 scale for the y-axis.
#' @param alpha Transparency level for the violin fill.
#' @param color Color for the violin outline.
#' @param draw_quantiles Numeric vector of quantiles to draw on the violin plot (default is 0.5 for median).
#' @param facet_var Optional variable to facet the plot by.
#' @param use_grid Logical indicating whether to use a grid layout for the plot.
#' @param use_jitter Logical indicating whether to add jittered data points.
#' @param jitter_size Size of the jittered points.
#' @param jitter_alpha Transparency level for the jittered points.
#'
#' @return A ggplot object representing the violin plot.
#'
#' @export
#'
plot_violin <- function(
  plot_data,
  x,
  y,
  fill = NULL,
  palette = NULL,
  title = NULL,
  subtitle = NULL,
  x_label = NULL,
  y_label = NULL,
  fill_label = NULL,
  summarize = TRUE,
  round = 1,
  use_pct = FALSE,
  use_1k = FALSE,
  expand = NULL,
  use_log10 = FALSE,
  alpha = 1,
  color = NA,
  draw_quantiles = 0.5,
  facet_var = NULL,
  use_grid = FALSE,
  use_jitter = FALSE,
  jitter_size = 0.1,
  jitter_alpha = 0.5
) {
  pixelatorR:::assert_class(plot_data, "data.frame")
  pixelatorR:::assert_single_value(x, type = "string", allow_null = TRUE)
  pixelatorR:::assert_single_value(y, type = "string", allow_null = TRUE)
  pixelatorR:::assert_single_value(fill, type = "string", allow_null = TRUE)
  pixelatorR:::assert_vector(palette, "character", allow_null = TRUE)
  pixelatorR:::assert_single_value(title, type = "string", allow_null = TRUE)
  pixelatorR:::assert_single_value(subtitle, type = "string", allow_null = TRUE)
  pixelatorR:::assert_single_value(x_label, type = "string", allow_null = TRUE)
  pixelatorR:::assert_single_value(y_label, type = "string", allow_null = TRUE)
  pixelatorR:::assert_single_value(fill_label, type = "string", allow_null = TRUE)
  pixelatorR:::assert_single_value(summarize, type = "bool")
  pixelatorR:::assert_single_value(round, type = "integer")
  pixelatorR:::assert_single_value(use_pct, type = "bool")
  pixelatorR:::assert_single_value(use_1k, type = "bool")
  pixelatorR:::assert_vector(expand, "numeric", n = 2, allow_null = TRUE)
  pixelatorR:::assert_single_value(use_log10, type = "bool")
  pixelatorR:::assert_single_value(alpha, type = "numeric")
  pixelatorR:::assert_within_limits(alpha, c(0, 1))
  pixelatorR:::assert_single_value(facet_var, type = "string", allow_null = TRUE)
  pixelatorR:::assert_single_value(use_grid, type = "bool")
  pixelatorR:::assert_single_value(use_jitter, type = "bool")
  pixelatorR:::assert_single_value(jitter_size, type = "numeric")
  pixelatorR:::assert_single_value(jitter_alpha, type = "numeric")
  pixelatorR:::assert_within_limits(jitter_alpha, c(0, 1))

  if (is.null(expand)) {
    expand <- waiver()
  } else {
    expand <- expansion(mult = expand)
  }

  p <- plot_data %>%
    {
      if (!is.null(fill)) {
        ggplot(., aes(!!sym(x), !!sym(y), fill = !!sym(fill))) +
          geom_violin(
            position = position_dodge(width = 0.9),
            alpha = alpha,
            color = color,
            scale = "width",
            drop = FALSE
          )
      } else {
        ggplot(., aes(!!sym(x), !!sym(y))) +
          geom_violin(
            position = position_dodge(width = 0.9),
            alpha = alpha,
            color = color,
            scale = "width",
            fill = "#DAD6D7",
            drop = FALSE
          )
      }
    }

  if (!is.null(facet_var)) {
    p <- p + facet_grid(as.formula(paste("~", facet_var)))
  }
  if (use_pct && !use_log10) {
    p <- p +
      scale_y_continuous(labels = scales::percent, expand = expand)
  }
  if (!use_pct && !use_log10) {
    p <- p +
      scale_y_continuous(expand = expand)
  }
  if (summarize) {
    plot_data_sum <-
      plot_data %>%
      {
        if (!is.null(fill)) {
          group_by(., !!sym(x), !!sym(fill))
        } else {
          group_by(., !!sym(x))
        }
      } %>%
      summarise(!!sym(paste0("median_", y)) := median(!!sym(y)), .groups = "drop")
    p <- p +
      geom_text(
        data = plot_data_sum,
        if (use_pct) {
          aes(
            y = !!sym(paste0("median_", y)),
            label = paste0(round(!!sym(paste0("median_", y)), round) * 100, "%")
          )
        } else if (use_1k) {
          aes(
            y = !!sym(paste0("median_", y)),
            label = paste0(round(!!sym(paste0("median_", y)) / 1000, round), "k")
          )
        } else {
          aes(
            y = !!sym(paste0("median_", y)),
            label = round(!!sym(paste0("median_", y)), round)
          )
        },
        position = position_dodge(width = 0.9),
        vjust = -0.8,
        size = 3
      )
  }
  if (!is.null(fill)) {
    p <- p +
      scale_fill_manual(values = palette)
  }
  if (!is.null(draw_quantiles)) {
    p <-
      p +
      draw_quantiles(p,
        draw_quantiles = draw_quantiles,
        facet_var = facet_var
      )
  }
  if (use_log10) {
    p <- p +
      scale_y_log10(expand = expand)
  }


  p <- p +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 60, hjust = 1),
      panel.grid = if (use_grid) element_line() else element_blank(),
    ) +
    labs(
      x = x_label,
      y = y_label,
      fill = fill_label,
      title = title,
      subtitle = subtitle
    )

  if (use_jitter) {
    p <- p +
      geom_jitter(
        size = jitter_size,
        position = position_jitter(height = 0),
        alpha = jitter_alpha
      )
  }

  return(p)
}



#' Draw quantiles on a violin plot
#'
#' This function adds quantile lines to a violin plot created with ggplot2.
#' It calculates the quantiles for each violin and draws horizontal lines at those quantiles.
#'
#' @param p A ggplot object representing the violin plot.
#' @param draw_quantiles A numeric vector of quantiles to draw (default is 0.5, which represents the median).
#' @param linewidth The width of the quantile lines (default is 1).
#' @param color The color of the quantile lines (default is "gray20").
#' @param facet_var An optional variable to facet the plot by (default is NULL).
#' @param use_real_quantile Logical indicating whether to use the actual quantile values (default is TRUE).
#' @param ... Additional arguments passed to `geom_segment`.
#'
#' @return A ggplot object with quantile lines added to the violin plot.
#'
#' @export
#'
draw_quantiles <-
  function(
    p,
    draw_quantiles = 0.5,
    linewidth = 1,
    color = "gray20",
    facet_var = NULL,
    use_real_quantile = TRUE,
    ...
  ) {
    pixelatorR:::assert_class(p, "ggplot")
    pixelatorR:::assert_single_value(draw_quantiles, "numeric", allow_null = TRUE)
    pixelatorR:::assert_single_value(linewidth, "numeric")
    pixelatorR:::assert_within_limits(linewidth, c(0, Inf))
    pixelatorR:::assert_single_value(color, "string")
    pixelatorR:::assert_valid_color(color)
    pixelatorR:::assert_single_value(facet_var, "string", allow_null = TRUE)
    pixelatorR:::assert_single_value(use_real_quantile, "bool")

    build <- ggplot_build(p)
    violin_data <- build$data[[1]]

    quantile_data <-
      violin_data %>%
      group_by(x, PANEL) %>%
      group_map(~ {
        .x <-
          bind_cols(.y, .x)
        # Get density distribution
        dens <- cumsum(.x$density) / sum(.x$density)
        # Create approximate cumulative density to actual density function
        ecdf <- stats::approxfun(dens, .x$y, ties = "ordered")

        # Get quantile
        if (use_real_quantile) {
          # Use the actual quantile values
          mask <-
            unclass(factor(p$data[[as_label(p$mapping$x)]])) == .y$x

          if (!is.null(facet_var)) {
            mask <-
              mask &
                (unclass(factor(p$data[[facet_var]])) == .y$PANEL)
          }

          ys <-
            stats::quantile(
              p$data[[as_label(p$mapping$y)]][mask],
              probs = draw_quantiles,
              na.rm = TRUE
            )
        } else {
          # Use the quantiles from the density distribution
          ys <- ecdf(draw_quantiles)
        }


        # Get violin width at quantile and coordinates
        bind_cols(
          .y,
          tibble(
            y = ys,
            violinwidth = with(
              .x,
              stats::approxfun(
                y,
                violinwidth *
                  (xmax - xmin)
              )
            )(ys)
          ) %>%
            mutate(
              xmin = .x$x[1] - violinwidth / 2,
              xmax = .x$x[1] + violinwidth / 2,
            )
        )
      }) %>%
      bind_rows() %>%
      select(x, y, PANEL, xmin, xmax)

    if (!is.null(facet_var)) {
      # If a facet variable is provided, add it to the quantile data
      quantile_data <-
        quantile_data %>%
        mutate(facet = factor(
          levels(p$data[[facet_var]])[PANEL],
          levels(p$data[[facet_var]])
        )) %>%
        rename(!!facet_var := facet)
    }

    geom_segment(
      data = quantile_data,
      aes(x = xmin, y = y, xend = xmax, yend = y),
      inherit.aes = FALSE,
      linewidth = linewidth,
      color = color,
      ...
    )
  }

#' Set sample levels in a data frame
#'
#' This function sets the levels of the `sample_alias` factor in a data frame.
#'
#' @param df A data frame containing a `sample_alias` column.
#' @param sample_levels A vector of sample levels to set for the `sample_alias` column (default is NULL).
#'
#' @return A data frame with the `sample_alias` column set to the specified levels.
#'
#' @export
#'
set_sample_levels <-
  function(
    df,
    sample_levels = NULL
  ) {
    pixelatorR:::assert_class(df, "data.frame")
    pixelatorR:::assert_vector(sample_levels, "character", n = 1, allow_null = TRUE)

    if (!is.null(sample_levels)) {
      df <-
        df %>%
        mutate(sample_alias = factor(sample_alias, sample_levels))
    } else {
      df <-
        df %>%
        mutate(sample_alias = factor(sample_alias))
    }

    return(df)
  }
