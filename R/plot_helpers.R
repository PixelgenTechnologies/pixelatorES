plot_anchor_slug <- function(...) {
  x <- paste(..., collapse = "-")
  x <- tolower(x)
  x <- gsub("[^a-z0-9]+", "-", x)
  gsub("^-+|-+$", "", x)
}

.tabset_state <- new.env(parent = emptyenv())
.tabset_state$depth <- 0L
.tabset_state$stack <- integer(0)

#' Reset tracked Quarto tabset depth
#'
#' Intended for tests. Report rendering resets state when
#' [register_tabset_chunk_hooks()] is called.
#'
#' @return `NULL`, invisibly.
#'
#' @noRd
reset_tabset_state <- function() {
  .tabset_state$depth <- 0L
  .tabset_state$stack <- integer(0)
  invisible(NULL)
}

#' Current number of opened tabsets that have not been closed
#'
#' @return Integer depth.
#'
#' @noRd
tabset_depth <- function() {
  .tabset_state$depth
}

#' Markdown fences that close leftover tabsets down to a depth
#'
#' Updates the tracked depth. Does not print; callers that run during knitr
#' chunk evaluation should `cat()` the result, while knitr output hooks should
#' return it so it is appended to the chunk markdown.
#'
#' @param until Integer depth to restore. Tabsets opened above this depth are
#'   closed.
#'
#' @return A string of `:::` fences, or `""` when nothing is unclosed.
#'
#' @noRd
unclosed_tabset_fences <- function(until = 0L) {
  until <- as.integer(until)
  if (length(until) != 1L || is.na(until) || until < 0L) {
    until <- 0L
  }

  n <- .tabset_state$depth - until
  if (n <= 0L) {
    return("")
  }

  .tabset_state$depth <- until
  paste(rep(":::", n), collapse = "\n")
}

#' Open a Quarto panel tabset
#'
#' Emits the fenced div that starts a `{.panel-tabset}` and records it so
#' leftover tabsets can be closed if a knitr chunk errors. Prefer this over
#' writing the fence with `cat()` in `.qmd` files.
#'
#' @return The new tabset depth, invisibly.
#'
#' @export
open_tabset <- function() {
  cat("::: {.panel-tabset .nav-pills}\n")
  .tabset_state$depth <- .tabset_state$depth + 1L
  invisible(.tabset_state$depth)
}

#' Close Tabset
#'
#' Closes a tabset opened with [open_tabset()] or [tabset_figure_table()].
#' Does nothing when no tracked tabset is open.
#'
#' @return The remaining tabset depth, invisibly.
#'
#' @export
close_tabset <- function() {
  fence <- unclosed_tabset_fences(until = max(0L, .tabset_state$depth - 1L))
  if (nzchar(fence)) {
    cat(fence, "\n", sep = "")
  }
  invisible(.tabset_state$depth)
}

#' Register knitr hooks that close leftover panel tabsets
#'
#' When Quarto is run with `execute: error: true`, a chunk can fail after
#' [open_tabset()] (or a `tabset_*` helper) has written an opening fence.
#' This hook restores the tabset depth from the start of the chunk by appending
#' the missing `:::` fences so later sections of the report still parse.
#'
#' Call once during report setup (see `inst/quarto/shared/preprocessing.qmd`).
#'
#' @return `NULL`, invisibly.
#'
#' @export
register_tabset_chunk_hooks <- function() {
  reset_tabset_state()
  knitr::opts_chunk$set(pixelatorES_tabset_guard = TRUE)
  knitr::knit_hooks$set(
    pixelatorES_tabset_guard = function(before, options, envir) {
      if (before) {
        .tabset_state$stack <- c(.tabset_state$stack, .tabset_state$depth)
        return(NULL)
      }

      n_stack <- length(.tabset_state$stack)
      start <- if (n_stack > 0L) {
        .tabset_state$stack[[n_stack]]
      } else {
        0L
      }
      if (n_stack > 0L) {
        .tabset_state$stack <- .tabset_state$stack[-n_stack]
      }

      fences <- unclosed_tabset_fences(until = start)
      if (!nzchar(fences)) {
        return(NULL)
      }

      paste0("\n\n", fences, "\n\n")
    }
  )
  invisible(NULL)
}

resolve_anchor_prefix <- function(anchor_prefix = NULL) {
  if (is.null(anchor_prefix)) {
    label <- knitr::opts_current$get("label")
    if (is.null(label) || !nzchar(label)) {
      return(NULL)
    }
    return(plot_anchor_slug(label))
  }
  anchor_prefix
}

.psock_pool <- new.env(parent = emptyenv())
.psock_pool$cl <- NULL
.psock_pool$n_workers <- 0L
.psock_pool$pkg_path <- NULL

.get_psock_cluster <- function(mc_cores) {
  mc_cores <- as.integer(mc_cores)
  if (is.na(mc_cores) || mc_cores <= 1L) {
    return(NULL)
  }

  ns_path <- tryCatch(getNamespaceInfo("pixelatorES", "path"), error = function(e) NULL)
  existing_cl <- .psock_pool$cl
  if (!is.null(existing_cl) && .psock_pool$n_workers == mc_cores && identical(.psock_pool$pkg_path, ns_path)) {
    alive <- tryCatch(
      {
        parallel::clusterCall(existing_cl, function() TRUE)
        TRUE
      },
      error = function(e) FALSE
    )
    if (alive) {
      return(existing_cl)
    }
  }

  if (!is.null(existing_cl)) {
    try(parallel::stopCluster(existing_cl), silent = TRUE)
    .psock_pool$cl <- NULL
    .psock_pool$n_workers <- 0L
  }

  cl <- tryCatch(
    {
      c_obj <- parallel::makePSOCKcluster(mc_cores)
      is_dev_source <- !is.null(ns_path) && file.exists(file.path(ns_path, "R", "plot_helpers.R"))
      parallel::clusterCall(
        c_obj,
        function(path, is_dev) {
          Sys.setenv(
            OMP_NUM_THREADS = "1",
            OPENBLAS_NUM_THREADS = "1",
            MKL_NUM_THREADS = "1",
            VECLIB_MAXIMUM_THREADS = "1",
            MAGICK_THREAD_LIMIT = "1"
          )
          if (is_dev && requireNamespace("pkgload", quietly = TRUE)) {
            suppressMessages(pkgload::load_all(path, quiet = TRUE))
          } else {
            suppressPackageStartupMessages(library(pixelatorES))
          }
          TRUE
        },
        ns_path,
        is_dev_source
      )
      c_obj
    },
    error = function(e) {
      NULL
    }
  )

  if (is.null(cl)) {
    return(NULL)
  }

  .psock_pool$cl <- cl
  .psock_pool$n_workers <- mc_cores
  .psock_pool$pkg_path <- ns_path
  cl
}

stop_psock_cluster <- function() {
  if (!is.null(.psock_pool$cl)) {
    try(parallel::stopCluster(.psock_pool$cl), silent = TRUE)
    .psock_pool$cl <- NULL
    .psock_pool$n_workers <- 0L
    .psock_pool$pkg_path <- NULL
  }
  invisible(NULL)
}

.build_single_abundance_marker_plot <- function(m, spec) {
  vals <- unname(spec$mat[m, ])
  marker_fac <- factor(m, levels = spec$ordered_markers)
  no_fill_guide <- ggplot2::guides(fill = "none")

  g_p1 <- spec$base_p1
  g_p1$marker <- marker_fac
  g_p1$normcount <- vals

  p1 <-
    plot_violin(
      g_p1,
      x = "sample_alias",
      y = "normcount",
      fill = "condition",
      fill_label = "Condition",
      title = m,
      y_label = "Normalized counts",
      summarize = FALSE,
      palette = spec$sample_palette,
      use_jitter = TRUE,
      use_grid = TRUE
    ) +
    no_fill_guide

  if (!spec$per_celltype) {
    p <- p1
  } else {
    g_p2 <- spec$base_p2
    g_p2$marker[!is.na(spec$p2_row_idx)] <- marker_fac
    g_p2$normcount <- vals[spec$p2_row_idx]

    p2 <-
      plot_violin(
        g_p2,
        x = "sample_alias",
        y = "normcount",
        fill = "condition",
        y_label = "Normalized counts",
        summarize = FALSE,
        palette = spec$sample_palette,
        facet_var = "celltype",
        use_jitter = TRUE,
        use_grid = TRUE
      ) +
      no_fill_guide

    p <- p1 / p2
  }

  if (!is.null(spec$positive_rate)) {
    pos_spec <- spec$positive_rate
    mrkr_data <- pos_spec$groups[[m]]
    if (!is.null(mrkr_data)) {
      p_pos <-
        mrkr_data |>
        ggplot2::ggplot(ggplot2::aes(sample_alias, celltype, fill = positive_rate)) +
        ggplot2::geom_tile(show.legend = FALSE) +
        ggplot2::geom_text(ggplot2::aes(label = round(positive_rate, 0)), size = 3) +
        ggplot2::scale_fill_gradientn(
          colours = pos_spec$palette,
          limits = c(0, 100)
        ) +
        ggplot2::scale_x_discrete(expand = ggplot2::expansion(0)) +
        ggplot2::scale_y_discrete(expand = ggplot2::expansion(0)) +
        ggplot2::theme_bw() +
        ggplot2::theme(
          axis.text.x = ggplot2::element_text(angle = 60, hjust = 1),
          axis.text.y = ggplot2::element_text(angle = 90, hjust = 0.5),
          panel.grid = ggplot2::element_blank(),
          axis.title.x = ggplot2::element_blank()
        ) +
        ggplot2::labs(
          title = paste0("Positive rate for ", unique(mrkr_data$marker)),
          subtitle = paste0("% Cells with CLR >= ", pos_spec$positive_threshold),
          y = "Cell type"
        )

      join_op <- pos_spec$join
      sizes <- pos_spec$sizes
      joined <- switch(
        join_op,
        "/" = p / p_pos,
        "|" = p | p_pos,
        "+" = p + p_pos
      )
      if (length(unique(sizes)) > 1L) {
        panel_sizes <- c(rep(sizes[1], length(joined) - 1L), sizes[2])
        if (join_op == "/") {
          joined <- joined + patchwork::plot_layout(heights = panel_sizes)
        } else if (join_op == "|") {
          joined <- joined + patchwork::plot_layout(widths = panel_sizes)
        }
      }
      p <- joined
    }
  }

  p
}

.build_single_proximity_marker_plot <- function(g_data, spec) {
  no_fill_guide <- ggplot2::guides(fill = "none")
  g_p1 <- g_data[spec$p1_idx, spec$p1_cols, drop = FALSE]
  g_p1$sample_alias <- spec$p1_sample_alias

  p1 <-
    plot_violin(
      g_p1,
      x = "sample_alias",
      y = "log2_ratio",
      fill = "condition",
      title = as.character(g_data$marker_1[1]),
      y_label = "Log2 ratio Proximity Score",
      summarize = FALSE,
      palette = spec$sample_palette,
      use_jitter = TRUE,
      use_grid = TRUE,
      hline = 0
    ) +
    no_fill_guide

  if (!spec$per_celltype) {
    return(p1)
  }

  g_p2 <- g_data[spec$p2_idx, spec$p2_cols, drop = FALSE]
  g_p2$sample_alias <- spec$p2_sample_alias
  g_p2$celltype <- spec$p2_celltype

  p2 <-
    plot_violin(
      g_p2,
      x = "sample_alias",
      y = "log2_ratio",
      fill = "condition",
      y_label = "Log2 ratio Proximity Score",
      summarize = FALSE,
      palette = spec$sample_palette,
      facet_var = "celltype",
      use_jitter = TRUE,
      use_grid = TRUE,
      hline = 0
    ) +
    no_fill_guide

  p1 / p2
}

.build_single_kit_qc_lot_plot <- function(g_data, spec) {
  blend_fn <- getExportedValue("ggblend", "blend")
  plot <- g_data |>
    ggplot2::ggplot(ggplot2::aes(x = condition, y = value, fill = lot_role)) +
    blend_fn(
      ggplot2::geom_violin(
        position = "identity",
        alpha = 1,
        color = NA,
        scale = "width"
      ),
      blend = "exclusion"
    ) +
    ggplot2::geom_violin(
      ggplot2::aes(group = interaction(condition, lot_role)),
      position = "identity",
      quantiles = 0.5,
      quantile.color = "black",
      quantile.linetype = 1L,
      color = NA,
      fill = NA,
      scale = "width"
    ) +
    ggplot2::scale_fill_manual(
      values = c(reference = "#DAD6D7", qc = "#C86584"),
      name = "Lot role"
    ) +
    ggplot2::labs(title = as.character(g_data$marker[1]), y = spec$y_label, x = NULL) +
    .violin_theme_nogrid

  if (!is.null(spec$hline)) {
    plot <- plot + ggplot2::geom_hline(yintercept = spec$hline, linetype = "dashed")
  }
  plot@plot_env <- emptyenv()
  plot
}

.render_worker_chunk <- function(chunk) {
  spec <- chunk$spec
  png_paths <- chunk$png_paths
  webp_paths <- chunk$webp_paths
  width <- chunk$width
  height <- chunk$height
  dpi <- chunk$dpi
  dev_args <- chunk$dev_args
  convert_webp <- chunk$convert_webp
  dev_mode <- chunk$dev_mode

  for (i in seq_along(png_paths)) {
    png_file <- png_paths[[i]]
    dir.create(dirname(png_file), recursive = TRUE, showWarnings = FALSE)
    args <- c(
      list(filename = png_file, width = width, height = height, units = "in", res = dpi),
      dev_args
    )
    do.call(grDevices::png, args)
    p <- try(
      if (identical(spec$type, "abundance_per_marker")) {
        .build_single_abundance_marker_plot(chunk$markers[[i]], spec)
      } else if (identical(spec$type, "kit_qc_lot_distribution")) {
        .build_single_kit_qc_lot_plot(chunk$groups[[i]], spec)
      } else {
        .build_single_proximity_marker_plot(chunk$groups[[i]], spec)
      },
      silent = TRUE
    )
    if (!inherits(p, "try-error")) {
      msg <- try(print(p), silent = TRUE)
    } else {
      msg <- p
    }
    if (inherits(msg, "try-error")) {
      if (dev_mode) {
        grDevices::dev.off()
        stop(msg)
      }
      print(ggplot2::ggplot() + ggplot2::geom_blank())
    }
    grDevices::dev.off()

    if (convert_webp) {
      img <- magick::image_read(png_file)
      magick::image_write(img, webp_paths[[i]], format = "webp", compression = "WebP")
      file.remove(png_file)
    }
  }
  if (convert_webp) webp_paths else png_paths
}

.parallel_render_plotlist <- function(spec, n_plots, mc_cores) {
  label <- knitr::opts_current$get("label")
  if (is.null(label) || !nzchar(label) || n_plots <= 1L) {
    return(NULL)
  }

  n_workers <- min(as.integer(mc_cores), n_plots)
  cl <- .get_psock_cluster(n_workers)
  if (is.null(cl)) {
    return(NULL)
  }

  opts <- knitr::opts_current$get()
  width <- if (!is.null(opts$fig.width)) opts$fig.width[1L] else 7
  height <- if (!is.null(opts$fig.height)) opts$fig.height[1L] else 7
  dpi <- if (!is.null(opts$dpi)) opts$dpi[1L] else 72
  dev_args <- opts$dev.args
  if (is.list(dev_args) && !is.null(dev_args$png)) {
    dev_args <- dev_args$png
  }
  if (!is.list(dev_args)) {
    dev_args <- list()
  }

  rel_png_paths <- vapply(
    seq_len(n_plots),
    function(i) as.character(knitr::fig_chunk(label, "png", i)),
    character(1)
  )
  rel_webp_paths <- sub("\\.png$", ".webp", rel_png_paths)

  base_dir <- knitr::opts_knit$get("base.dir")
  abs_png_paths <- if (is.character(base_dir) && nzchar(base_dir)) {
    file.path(base_dir, rel_png_paths)
  } else {
    normalizePath(rel_png_paths, mustWork = FALSE)
  }
  abs_webp_paths <- sub("\\.png$", ".webp", abs_png_paths)

  plot_hook <- knitr::knit_hooks$get("plot")
  convert_webp <- identical(plot_hook, convert_png_to_webp)
  dev_mode <- isTRUE(getOption("pixelatorES.dev_mode", default = FALSE))

  idx_splits <- parallel::splitIndices(n_plots, n_workers)
  idx_splits <- idx_splits[lengths(idx_splits) > 0L]

  chunks <- lapply(idx_splits, function(idx) {
    if (identical(spec$type, "abundance_per_marker")) {
      sub_markers <- spec$ordered_markers[idx]
      sub_spec <- spec
      sub_spec$mat <- spec$mat[sub_markers, , drop = FALSE]
      if (!is.null(spec$positive_rate)) {
        sub_spec$positive_rate <- spec$positive_rate
        sub_spec$positive_rate$groups <- spec$positive_rate$groups[sub_markers]
      }
      list(
        spec = sub_spec,
        markers = sub_markers,
        png_paths = abs_png_paths[idx],
        webp_paths = abs_webp_paths[idx],
        width = width,
        height = height,
        dpi = dpi,
        dev_args = dev_args,
        convert_webp = convert_webp,
        dev_mode = dev_mode
      )
    } else {
      sub_spec <- spec
      sub_spec$groups <- NULL
      list(
        spec = sub_spec,
        groups = spec$groups[idx],
        png_paths = abs_png_paths[idx],
        webp_paths = abs_webp_paths[idx],
        width = width,
        height = height,
        dpi = dpi,
        dev_args = dev_args,
        convert_webp = convert_webp,
        dev_mode = dev_mode
      )
    }
  })

  parallel::parLapply(cl, chunks, function(ch) {
    pixelatorES:::.render_worker_chunk(ch)
  })

  if (convert_webp) rel_webp_paths else rel_png_paths
}

#' Header wrap plots
#'
#' Wraps a list of plots with headers for each plot.
#'
#' @param plots List of plots to be wrapped.
#' @param level Integer indicating the header level for the headers (default is 2).
#' @param anchor_prefix Prefix for in-report deep-link anchor IDs. `NULL` (default)
#'   auto-detects from the current knitr chunk label; pass a string to override.
#' @param mc_cores Optional integer number of cores for parallel plot rendering
#'   during knitr/Quarto execution. When `NULL`, uses the `mc_cores` attached to
#'   `plots` (from `es_data$params$mc_cores`) or `getOption("pixelatorES.mc_cores", 1L)`.
#'
#' @return A formatted list of plots with headers.
#'
#' @export
title_plotlist <- function(plots, level = 2, anchor_prefix = NULL, mc_cores = NULL) {
  pixelatorR:::assert_class(plots, "list")
  for (plot in plots) pixelatorR:::assert_class(plot, c("ggplot", "datatables"))

  prefix <- resolve_anchor_prefix(anchor_prefix)
  nams <- names(plots)

  if (is.null(nams)) nams <- rep(" ", length(plots))

  spec <- attr(plots, "parallel_render_spec")
  if (is.null(mc_cores)) {
    mc_cores <- if (!is.null(spec$mc_cores)) {
      spec$mc_cores
    } else {
      getOption("pixelatorES.mc_cores", default = 1L)
    }
  }
  mc_cores <- as.integer(mc_cores)

  rel_paths <- NULL
  if (!is.null(spec) && !is.na(mc_cores) && mc_cores > 1L && isTRUE(getOption("knitr.in.progress"))) {
    rel_paths <- .parallel_render_plotlist(spec, length(plots), mc_cores)
  }

  if (!is.null(rel_paths)) {
    opts <- knitr::opts_current$get()
    plot_hook <- knitr::knit_hooks$get("plot")
    use_md_hook <- identical(plot_hook, convert_png_to_webp) || is.null(plot_hook)
    for (tab in seq_along(plots)) {
      cat(paste0(strrep("#", level), " ", nams[tab], "\n\n"))

      if (!is.null(prefix)) {
        cat(sprintf(
          '<div class="plot-anchor" data-anchor-id="%s" aria-hidden="true"></div>\n\n',
          plot_anchor_slug(prefix, nams[tab])
        ))
      }

      md_str <- if (use_md_hook) {
        knitr::hook_plot_md(rel_paths[[tab]], opts)
      } else {
        plot_hook(rel_paths[[tab]], opts)
      }
      cat(md_str)
      if (!grepl("\n$", md_str)) {
        cat("\n\n")
      }
    }
    return(invisible(NULL))
  }

  # Loop through each plot and create a tab for each
  for (tab in seq_along(plots)) {
    # Generate a header for each plot in the tabset
    cat(paste0(strrep("#", level), " ", nams[tab], "\n\n"))

    if (!is.null(prefix)) {
      cat(sprintf(
        '<div class="plot-anchor" data-anchor-id="%s" aria-hidden="true"></div>\n\n',
        plot_anchor_slug(prefix, nams[tab])
      ))
    }

    # Print the plot
    msg <- try(print(plots[[tab]]), silent = TRUE)
    if (inherits(msg, "try-error")) {
      if (getOption("pixelatorES.dev_mode", default = FALSE)) stop(msg)

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
#' @param close Logical indicating whether to close the tabset after plotting (default is TRUE).
#' @param anchor_prefix Prefix for in-report deep-link anchor IDs. `NULL` (default)
#'   auto-detects from the current knitr chunk label; pass a string to override.
#' @param mc_cores Optional integer number of cores for parallel plot rendering
#'   during knitr/Quarto execution.
#'
#' @return A formatted tabset of plots.
#'
#' @export
tabset_plotlist <-
  function(
    plots,
    level = 2,
    close = TRUE,
    anchor_prefix = NULL,
    mc_cores = NULL
  ) {
    pixelatorR:::assert_class(plots, "list")
    for (plot in plots) pixelatorR:::assert_class(plot, c("ggplot", "datatables"))

    prefix <- resolve_anchor_prefix(anchor_prefix)

    open_tabset()
    if (close) {
      on.exit(close_tabset(), add = TRUE)
    }

    title_plotlist(plots, level, anchor_prefix = prefix, mc_cores = mc_cores)
  }


#' Tabset nested plots
#'
#' Tabset a list of plots, which may contain nested lists of plots, at a given header level.
#'
#' @param plots List of plots or list of lists of plots to be tabsetted.
#' @param level Integer indicating the header level for the tabset title (default is 2).
#' @param close Logical indicating whether to close the tabset after plotting (default is TRUE).
#' @param anchor_prefix Prefix for in-report deep-link anchor IDs. `NULL` (default)
#'   auto-detects from the current knitr chunk label; pass a string to override.
#'
#' @return A formatted tabset of plots.
#'
#' @export
tabset_nested_plotlist <-
  function(
    plots,
    level = 2,
    close = TRUE,
    anchor_prefix = NULL
  ) {
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

    prefix <- resolve_anchor_prefix(anchor_prefix)

    cat("\n\n")
    open_tabset()
    if (close) {
      on.exit(close_tabset(), add = TRUE)
    }

    for (tab in seq_along(plots)) {
      if (inherits(plots[[tab]], "list")) {
        # If the plot is a list, call tabset_nested_plotlist recursively

        cat(paste0(strrep("#", level), " ", names(plots)[tab], "\n"))

        child_prefix <- if (!is.null(prefix)) {
          plot_anchor_slug(prefix, names(plots)[tab])
        } else {
          NULL
        }

        tabset_nested_plotlist(
          plots[[tab]],
          level + 1,
          anchor_prefix = child_prefix
        )
      } else if (inherits(plots[[tab]], "ggplot")) {
        # Otherwise, just print the plot
        title_plotlist(
          list(plots[[tab]]) %>% set_names(names(plots)[tab]),
          level,
          anchor_prefix = prefix
        )
      }
    }
  }


#' Tabset a figure and a table
#'
#' Tabset a figure and a table at a given header level. After using this
#' function, you should close the tabset with [close_tabset()]. If the chunk
#' errors before that call, [register_tabset_chunk_hooks()] closes the fence.
#'
#' @param figure A ggplot object representing the figure to be tabsetted or a list of ggplot objects to be tabsetted as
#' nested tabs.
#' @param table A table to be tabsetted alongside the figure.
#' @param level Integer indicating the header level for the tabset title (default is 2).
#' @param mode Character indicating the mode of the tabset ("tabset", "tabset_nested", or "title").
#'
#' @return A formatted tabset containing the figure and table.
#'
#' @export
tabset_figure_table <- function(figure, table, level = 2,
                                mode = "tabset") {
  pixelatorR:::assert_class(table, "datatables")
  mode <- match.arg(mode, c("tabset", "tabset_nested", "title"))

  set_func <-
    switch(mode,
      tabset = tabset_plotlist,
      tabset_nested = tabset_nested_plotlist,
      title = title_plotlist
    )

  open_tabset()

  cat(paste0(strrep("#", level), " Figure\n\n"))

  if (inherits(figure, "list")) {
    set_func(figure, level + 1)
  } else {
    pixelatorR:::assert_class(figure, "ggplot")
    print(figure)
  }
  cat("\n\n")
  cat(paste0(strrep("#", level), " Table\n\n"))

  return(table)
}

#' Create a Void Plot
#'
#' Creates a ggplot object with no data and a void theme.
#'
#' @return A ggplot object with a void theme.
#'
#' @export
#'
plot_void <-
  function() {
    ggplot() +
      theme_void() +
      theme(plot.margin = margin(0, 0, 0, 0))
  }

#' Extract Legend from a ggplot
#'
#' Extracts the legend from a ggplot object and returns it as a separate ggplot object.
#'
#' @param plot A ggplot object from which to extract the legend.
#' @param legend_position The position of the legend to extract (default is "bottom").
#' @param plot_height The height of the plot in inches (default is 5).
#' @param justify_top Logical indicating whether to justify the legend to the top (default is TRUE).
#'
#' @return A ggplot object containing only the legend.
#'
#' @export
#'
.extract_legend <-
  function(
    plot,
    legend_position = "bottom",
    plot_height = 5,
    justify_top = TRUE
  ) {
    pixelatorR:::assert_class(plot, "ggplot")
    pixelatorR:::assert_single_value(legend_position, type = "string")
    pixelatorR:::assert_single_value(plot_height, type = "numeric")
    pixelatorR:::assert_single_value(justify_top, type = "bool")
    legend_position <- match.arg(legend_position, c("bottom", "top", "left", "right", "none"))

    if (legend_position == "none") {
      return(plot_void())
    }

    plot <-
      plot +
      theme(
        legend.background = element_blank(),
        legend.box.background = element_blank(),
        legend.key = element_blank()
      )

    # Extract the legend from the plot
    g <- ggplotGrob(plot)

    guide_grob <- g$grobs[[which(g$layout$name == paste0("guide-box-", legend_position))]]

    legend_plot <-
      plot_void() +
      annotation_custom(guide_grob)

    if (justify_top) {
      # get legend height in mm (sum of the gtable row heights)
      legend_height <-
        sum(grid::convertUnit(guide_grob$heights, "in", valueOnly = TRUE))

      # create a column for the legend where the top row is the legend and the bottom row is spacer
      # use heights proportional to the absolute heights (so the legend gets the space it needs)
      legend_plot <-
        legend_plot /
        plot_spacer() +
        plot_layout(heights = c(legend_height, max(0, plot_height - legend_height)))
    }

    return(legend_plot)
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
#' @param extract_legend Logical indicating whether to extract the legend from the plot (default is FALSE). If TRUE,
#' the function will return a list containing the plot and the legend separately.
#' @param plot_height The height of the plot in inches (default is 5).
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
           legend_cols = 2,
           extract_legend = FALSE,
           plot_height = 5) {
    pixelatorR:::assert_class(object, "Seurat")
    pixelatorR:::assert_single_value(plot_reduction, type = "string")
    pixelatorR:::assert_vector(dims, "numeric", n = 2)
    pixelatorR:::assert_vector(metavars, "character", n = 1)
    pixelatorR:::assert_vector(pal, "character", n = 1, allow_null = TRUE)
    pixelatorR:::assert_single_value(label, type = "bool")
    pixelatorR:::assert_single_value(plot_title, type = "string")
    pixelatorR:::assert_single_value(xaxis_title, type = "bool")
    pixelatorR:::assert_single_value(yaxis_title, type = "bool")
    pixelatorR:::assert_single_value(fix_aspect_ratio, type = "bool")
    pixelatorR:::assert_single_value(legend_position, type = "string")
    pixelatorR:::assert_single_value(legend_cols, type = "numeric")
    pixelatorR:::assert_single_value(extract_legend, type = "bool")

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
        legend.key.size = unit(0.2, "cm"),
        legend.spacing.x = unit(0.2, "cm"),
        legend.direction = "vertical",
        legend.text = element_text(size = 8)
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

    if (extract_legend) {
      legend <- .extract_legend(p,
        legend_position = legend_position,
        plot_height = plot_height
      )
      p <- p + theme(legend.position = "none")

      return(list(plot = p, legend = legend))
    } else {
      return(p)
    }
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
#' @param fix_aspect_ratio Logical indicating whether to fix the aspect ratio of the plot (default is TRUE).
#' @param legend_position The position of the legend (default is "none").
#' @param legend_cols The number of columns in the legend (default is 2).
#'
#' @return A list of ggplot objects, each representing the embedding for a sample.
#' @export
plot_embeddings_samplewise <-
  function(
    samples,
    object,
    plot_reduction = "pca",
    dims = 1:2,
    metavars = "seurat_clusters",
    pal = NULL,
    label = TRUE,
    plot_title = "",
    fix_aspect_ratio = TRUE,
    legend_position = "none",
    legend_cols = 2
  ) {
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


    if (fix_aspect_ratio) {
      V1_range <- range(Embeddings(object, plot_reduction)[, dims[1]])
      V2_range <- range(Embeddings(object, plot_reduction)[, dims[2]])

      range_diff <- diff(V2_range) - diff(V1_range)

      if (range_diff < 0) {
        V2_range <- V2_range + c(-1, 1) * abs(range_diff) / 2
      } else if (range_diff > 0) {
        V1_range <- V1_range + c(-1, 1) * abs(range_diff) / 2
      }
    } else {
      V1_range <- NULL
      V2_range <- NULL
    }


    samples %>%
      set_names() %>%
      lapply(function(samp_id) {
        p <-
          object %>%
          subset(subset = sample_alias == samp_id) %>%
          plot_embedding(
            plot_reduction = plot_reduction,
            dims = dims,
            metavars = metavars,
            pal = pal,
            label = label,
            plot_title = plot_title,
            fix_aspect_ratio = FALSE,
            legend_position = legend_position,
            legend_cols = legend_cols
          )
        if (fix_aspect_ratio) {
          p <-
            p +
            scale_x_continuous(limits = V1_range) +
            scale_y_continuous(limits = V2_range)
        }
        p
      })
  }


# Pre-cached base violin themes to avoid repeated S7 theme construction/validation in ggplot2 >= 4.0
.violin_theme_grid <-
  ggplot2::theme_bw() +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 60, hjust = 1),
    panel.grid = ggplot2::element_line()
  )

.violin_theme_nogrid <-
  ggplot2::theme_bw() +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 60, hjust = 1),
    panel.grid = ggplot2::element_blank()
  )

# Vectorized 32-bit bit-reversal for base-2 Van der Corput sequence.
# Produces identical values to vipor:::vanDerCorput(n, base = 2, start = 1)
# without scalar R loops.
fast_vdc_base2 <- function(n, base = 2L, start = 1L) {
  if (n <= 0L) {
    return(numeric(0))
  }
  if (base != 2L) {
    return(getFromNamespace("vanDerCorput", "vipor")(n, base, start))
  }
  x <- seq.int(as.integer(start), length.out = as.integer(n))
  x <- bitwOr(bitwShiftL(bitwAnd(x, 0x55555555L), 1L), bitwAnd(bitwShiftR(x, 1L), 0x55555555L))
  x <- bitwOr(bitwShiftL(bitwAnd(x, 0x33333333L), 2L), bitwAnd(bitwShiftR(x, 2L), 0x33333333L))
  x <- bitwOr(bitwShiftL(bitwAnd(x, 0x0F0F0F0FL), 4L), bitwAnd(bitwShiftR(x, 4L), 0x0F0F0F0FL))
  x <- bitwOr(bitwShiftL(bitwAnd(x, 0x00FF00FFL), 8L), bitwAnd(bitwShiftR(x, 8L), 0x00FF00FFL))
  low16 <- bitwAnd(x, 0xFFFFL)
  high16 <- bitwAnd(bitwShiftR(x, 16L), 0xFFFFL)
  (as.numeric(low16) * 65536 + as.numeric(high16)) / 4294967296
}

fast_offset_single_group <- function(
  y,
  maxLength = NULL,
  method = "quasirandom",
  nbins = NULL,
  adjust = 1
) {
  if (!identical(method, "quasirandom")) {
    return(vipor::offsetSingleGroup(
      y,
      maxLength = maxLength,
      method = method,
      nbins = nbins,
      adjust = adjust
    ))
  }
  n <- length(y)
  if (n == 0L) {
    return(NULL)
  }
  if (n == 1L) {
    return(0)
  }
  if (is.null(nbins)) {
    nbins <- 1024L
  }
  subgroup_width <- if (is.null(maxLength) || maxLength <= 0) 1 else sqrt(n / maxLength)
  dens <- stats::density(y, n = nbins, adjust = adjust)
  dens$y <- dens$y / max(dens$y)
  offset <- fast_vdc_base2(n)[rank(y, ties.method = "first")]
  point_densities <- stats::approx(dens$x, dens$y, y)$y
  (offset - 0.5) * 2 * point_densities * subgroup_width
}

fast_offset_quasirandom <- function(
  data,
  width = 0.4,
  vary.width = FALSE,
  max.length = NULL,
  ...
) {
  x.offset <- vipor::aveWithArgs(
    data$y,
    data$x,
    FUN = fast_offset_single_group,
    maxLength = if (vary.width) max.length else NULL,
    ...
  )
  data$x <- data$x + x.offset * width
  if ("xend" %in% colnames(data) && "yend" %in% colnames(data)) {
    x.offset <- vipor::aveWithArgs(
      data$yend,
      data$xend,
      FUN = fast_offset_single_group,
      maxLength = if (vary.width) max.length else NULL,
      ...
    )
    data$xend <- data$xend + x.offset * width
  }
  data
}

PositionQuasirandomFast <- ggplot2::ggproto(
  "PositionQuasirandomFast",
  ggbeeswarm:::PositionQuasirandom,
  compute_panel = function(data, params, scales) {
    data <- ggplot2::flip_data(data, params$flipped_aes)
    data <- ggplot2:::collide(
      data,
      params$dodge.width,
      name = "position_quasirandom",
      strategy = ggplot2:::pos_dodge,
      check.width = FALSE
    )
    if (is.null(params$width)) {
      params$width <- ggplot2::resolution(data$x, zero = FALSE) * 0.4
    }
    if (!is.null(params$dodge.width)) {
      data <- split(data, data$group)
    } else {
      data <- split(data, data$x)
    }
    data <- lapply(
      data,
      fast_offset_quasirandom,
      method = params$method,
      width = params$width,
      vary.width = params$varwidth,
      adjust = params$bandwidth,
      nbins = params$nbins,
      max.length = params$max.length
    )
    data <- do.call(rbind, data)
    ggplot2::flip_data(data, params$flipped_aes)
  }
)

.position_qr_default <- ggplot2::ggproto(
  NULL,
  PositionQuasirandomFast,
  width = NULL,
  varwidth = FALSE,
  bandwidth = 0.5,
  nbins = NULL,
  method = "quasirandom",
  dodge.width = 0.9,
  na.rm = FALSE,
  orientation = NULL
)

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
#' @param use_jitter Logical indicating whether to add quasirandom data points.
#' @param jitter_size Size of the quasirandom points.
#' @param jitter_alpha Transparency level for the quasirandom points (default 0.35).
#' @param hline Optional horizontal line to draw on the plot.
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
  jitter_alpha = 0.35,
  hline = NULL
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
  pixelatorR:::assert_single_value(hline, type = "numeric", allow_null = TRUE)

  if (is.null(expand)) {
    expand <- waiver()
  } else {
    expand <- expansion(mult = expand)
  }

  # Check that there are enough datapoints to plot violins using fast integer tabulation
  y_vals <- plot_data[[y]]
  valid <- !is.na(y_vals)
  grp_cols <- unique(c(x, facet_var))
  if (length(grp_cols) == 1L) {
    g <- as.integer(as.factor(plot_data[[grp_cols]][valid]))
  } else {
    g <- as.integer(interaction(plot_data[valid, grp_cols, drop = FALSE], drop = TRUE))
  }
  counts <- tabulate(g)
  max_complete_datapoints <- if (length(counts) > 0L) max(counts) else 0L

  has_q <- !is.null(draw_quantiles) && length(draw_quantiles) > 0L
  has_all_facet_levels <- FALSE

  if (max_complete_datapoints <= 1L) {
    p <- ggplot(plot_data, aes(!!sym(x), !!sym(y))) +
      geom_hline(yintercept = hline)
  } else {
    q_cols <- unique(c(x, fill, facet_var))
    g_all <- if (length(q_cols) == 1L) {
      as.integer(as.factor(plot_data[[q_cols]]))
    } else {
      as.integer(interaction(plot_data[, q_cols, drop = FALSE], drop = TRUE))
    }
    q_counts <- tabulate(g_all[valid])
    has_all_facet_levels <- !is.null(facet_var) &&
      all(displayed_cell_types %in% unique(as.character(plot_data[[facet_var]])))

    violin_data <- if (has_all_facet_levels) {
      keep <- valid & !is.na(g_all) & q_counts[g_all] >= 2L
      if (all(keep)) NULL else plot_data[keep, , drop = FALSE]
    } else if (any(q_counts == 1L)) {
      keep <- !valid | (!is.na(g_all) & q_counts[g_all] >= 2L)
      plot_data[keep, , drop = FALSE]
    } else {
      NULL
    }

    if (!is.null(fill)) {
      p <- ggplot(plot_data, aes(!!sym(x), !!sym(y), fill = !!sym(fill))) +
        geom_hline(yintercept = hline) +
        geom_violin(
          data = violin_data,
          position = position_dodge(width = 0.9),
          alpha = alpha,
          color = color,
          scale = "width",
          drop = FALSE,
          quantiles = if (has_q) draw_quantiles else c(0.25, 0.5, 0.75),
          quantile.linetype = if (has_q) 1L else 0L,
          quantile.color = "black"
        )
    } else {
      p <- ggplot(plot_data, aes(!!sym(x), !!sym(y))) +
        geom_hline(yintercept = hline) +
        geom_violin(
          data = violin_data,
          position = position_dodge(width = 0.9),
          alpha = alpha,
          color = color,
          scale = "width",
          fill = "#DAD6D7",
          drop = FALSE,
          quantiles = if (has_q) draw_quantiles else c(0.25, 0.5, 0.75),
          quantile.linetype = if (has_q) 1L else 0L,
          quantile.color = "black"
        )
    }
    if (has_all_facet_levels) {
      p$layers[[length(p$layers)]]$geom_params$width <- 0.9
    }
  }

  if (!is.null(facet_var)) {
    p <- p + facet_grid(as.formula(paste("~", facet_var)))
  }
  if (has_all_facet_levels) {
    x_limits <- if (is.factor(plot_data[[x]])) {
      levels(plot_data[[x]])
    } else {
      unique(as.character(plot_data[[x]]))
    }
    p <- p + ggplot2::scale_x_discrete(limits = x_limits, drop = FALSE)
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
  if (use_log10) {
    p <- p +
      scale_y_log10(expand = expand)
  }

  p <- p +
    (if (use_grid) .violin_theme_grid else .violin_theme_nogrid) +
    labs(
      x = x_label,
      y = y_label,
      fill = fill_label,
      title = title,
      subtitle = subtitle
    )

  if (use_jitter) {
    point_group_aes <-
      if (!is.null(fill)) {
        aes(group = interaction(!!sym(x), !!sym(fill)))
      } else {
        aes(group = !!sym(x))
      }

    p <- p +
      ggplot2::layer(
        data = if (has_all_facet_levels && !all(valid)) plot_data[valid, , drop = FALSE] else NULL,
        mapping = point_group_aes,
        stat = "identity",
        geom = ggplot2::GeomPoint,
        position = .position_qr_default,
        show.legend = NA,
        inherit.aes = TRUE,
        params = list(
          na.rm = FALSE,
          color = "grey20",
          alpha = jitter_alpha,
          size = jitter_size
        )
      )
  }

  p@plot_env <- emptyenv()
  return(p)
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
      sample_levels <- sample_levels[sample_levels %in% unique(df$sample_alias)]

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
