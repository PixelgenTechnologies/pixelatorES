test_that("Tab setting and title setting work as expected", {
  pixelatorES:::reset_tabset_state()
  on.exit(pixelatorES:::reset_tabset_state(), add = TRUE)

  g <-
    ggplot() +
    theme_void()

  plot_list <-
    list("a" = g, "b" = g, "c" = g)

  expect_equal(
    capture.output(title_plotlist(plot_list, anchor_prefix = NULL)),
    c("## a", "", "", "", "## b", "", "", "", "## c", "", "", "")
  )
  expect_equal(
    capture.output(title_plotlist(plot_list, level = 3, anchor_prefix = NULL)),
    c("### a", "", "", "", "### b", "", "", "", "### c", "", "", "")
  )


  expect_equal(
    capture.output(tabset_plotlist(plot_list, anchor_prefix = NULL)),
    c(
      "::: {.panel-tabset .nav-pills}", "## a", "", "", "", "## b",
      "", "", "", "## c", "", "", "", ":::"
    )
  )
  expect_equal(
    capture.output(tabset_plotlist(plot_list, level = 3, anchor_prefix = NULL)),
    c(
      "::: {.panel-tabset .nav-pills}", "### a", "", "", "", "### b",
      "", "", "", "### c", "", "", "", ":::"
    )
  )

  # Capture plots that will throw an error

  ## Throw no error outside dev mode
  options(pixelatorES.dev_mode = FALSE)

  g_bad <-
    ggplot(tibble(), aes(notacolumn)) +
    theme_void()

  plot_list <-
    list("a" = g, "b" = g, "c" = g_bad)

  expect_no_error(title_plotlist(plot_list, anchor_prefix = NULL))

  ## Throw error in dev mode
  options(pixelatorES.dev_mode = TRUE)

  expect_error(title_plotlist(plot_list, anchor_prefix = NULL))
  expect_error(tabset_plotlist(plot_list, anchor_prefix = NULL))
  expect_equal(pixelatorES:::tabset_depth(), 0L)

  options(pixelatorES.dev_mode = FALSE)

  # Nested plot list

  plot_list <-
    list("a" = g, "b" = g, "c" = g)

  nested_plot_list <-
    list(
      "a" = g, "b" = g, "c" = g,
      plot_list = plot_list
    )

  expect_equal(
    capture.output(tabset_nested_plotlist(nested_plot_list, level = 3, anchor_prefix = NULL)),
    c(
      "", "", "::: {.panel-tabset .nav-pills}", "### a", "", "",
      "", "### b", "", "", "", "### c", "", "", "", "### plot_list",
      "", "", "::: {.panel-tabset .nav-pills}", "#### a", "", "", "",
      "#### b", "", "", "", "#### c", "", "", "", ":::", ":::"
    )
  )
})

test_that("plot anchor slugs are stable and URL-safe", {
  expect_equal(pixelatorES:::plot_anchor_slug("HLA-DR"), "hla-dr")
  expect_equal(
    pixelatorES:::plot_anchor_slug("abundance_per_marker", "CD3"),
    "abundance-per-marker-cd3"
  )
})

test_that("title_plotlist emits anchor divs when anchor_prefix is set", {
  g <- ggplot() + theme_void()
  plot_list <- list("CD3" = g, "HLA-DR" = g)

  out <- capture.output(
    title_plotlist(plot_list, level = 5, anchor_prefix = "abundance-per-marker")
  )

  expect_true(any(grepl('data-anchor-id="abundance-per-marker-cd3"', out)))
  expect_true(any(grepl('data-anchor-id="abundance-per-marker-hla-dr"', out)))
})

test_that("tabset_nested_plotlist extends anchor prefixes for nested tabs", {
  pixelatorES:::reset_tabset_state()
  on.exit(pixelatorES:::reset_tabset_state(), add = TRUE)

  g <- ggplot() + theme_void()
  nested_plot_list <- list(
    "B cell" = list("CD19" = g)
  )

  out <- capture.output(
    tabset_nested_plotlist(
      nested_plot_list,
      level = 3,
      anchor_prefix = "coloc-celltype"
    )
  )

  expect_true(any(grepl('data-anchor-id="coloc-celltype-b-cell-cd19"', out)))
})

test_that("tabset depth is tracked and leftover fences can be flushed", {
  pixelatorES:::reset_tabset_state()
  on.exit(pixelatorES:::reset_tabset_state(), add = TRUE)

  expect_equal(pixelatorES:::tabset_depth(), 0L)
  expect_equal(
    capture.output(open_tabset()),
    "::: {.panel-tabset .nav-pills}"
  )
  expect_equal(pixelatorES:::tabset_depth(), 1L)
  expect_equal(capture.output(close_tabset()), ":::")
  expect_equal(pixelatorES:::tabset_depth(), 0L)

  capture.output(open_tabset())
  capture.output(open_tabset())
  expect_equal(pixelatorES:::tabset_depth(), 2L)
  expect_equal(
    pixelatorES:::unclosed_tabset_fences(until = 0L),
    ":::\n:::"
  )
  expect_equal(pixelatorES:::tabset_depth(), 0L)
  expect_equal(capture.output(close_tabset()), character())
})

test_that("tabset chunk hook closes fences left open when a chunk errors", {
  pixelatorES:::reset_tabset_state()
  old_opt <- knitr::opts_chunk$get("pixelatorES_tabset_guard")
  old_hook <- knitr::knit_hooks$get("pixelatorES_tabset_guard")
  on.exit(
    {
      pixelatorES:::reset_tabset_state()
      knitr::opts_chunk$set(pixelatorES_tabset_guard = old_opt)
      if (!is.null(old_hook)) {
        knitr::knit_hooks$set(pixelatorES_tabset_guard = old_hook)
      }
    },
    add = TRUE
  )

  rmd <- tempfile(fileext = ".Rmd")
  on.exit(unlink(rmd), add = TRUE)

  writeLines(
    c(
      "```{r, results='asis', error=TRUE}",
      "open_tabset()",
      "cat('## Inside\\n\\n')",
      "stop('chunk failed after opening a tabset')",
      "close_tabset()",
      "```",
      "",
      "after the failed chunk"
    ),
    rmd
  )

  register_tabset_chunk_hooks()
  md_file <- knitr::knit(rmd, quiet = TRUE, envir = new.env(parent = globalenv()))
  on.exit(unlink(md_file), add = TRUE)
  md <- paste(readLines(md_file), collapse = "\n")

  expect_match(md, "::: \\{\\.panel-tabset \\.nav-pills\\}", perl = TRUE)
  expect_match(md, "after the failed chunk", perl = TRUE)
  open_at <- regexpr("::: \\{\\.panel-tabset", md)
  after_at <- regexpr("after the failed chunk", md)
  closings <- gregexpr("(^|\\n):::(\\n|$)", md, perl = TRUE)[[1]]
  expect_true(any(closings > open_at & closings < after_at))
  expect_equal(pixelatorES:::tabset_depth(), 0L)
})

test_that("Embedding plots work as expected", {
  set.seed(37)
  seur <-
    CreateSeuratObject(
      counts = matrix(
        c(
          rpois(1000, 40),
          rpois(1000, 5)
        )[sample(1:2000, 2000)],
        nrow = 100,
        ncol = 20,
        dimnames = list(
          paste0("Feature", 1:100),
          paste0("Cell", 1:20)
        )
      ) %>%
        as("dgCMatrix")
    ) %>%
    AddMetaData(
      metadata = data.frame(
        sample_alias = rep(c("A", "B"), each = 10),
        sample_type = rep(
          c("Unstimulated", "Stimulated"),
          each = 5,
          times = 2
        ),
        row.names = paste0("Cell", 1:20)
      )
    )

  seur <-
    seur %>%
    NormalizeData(verbose = FALSE) %>%
    ScaleData(verbose = FALSE) %>%
    RunPCA(npcs = 2, features = rownames(seur), verbose = FALSE) %>%
    AddMetaData(
      metadata = data.frame(
        seurat_clusters =
          sample(1:10, ncol(seur), replace = TRUE),
        row.names = paste0("Cell", 1:20)
      )
    )


  expect_no_error(p <- plot_embedding(seur))
  expect_s3_class(p, "ggplot")
  expect_no_error(p <- plot_embeddings_samplewise(c("A", "B"), seur))
  expect_type(p, "list")
  expect_equal(
    capture.output(p),
    c("$A", "", "$B", "")
  )
})

test_that("`plot_violin` works as expected", {
  plot_data <-
    structure(
      list(sample_alias = c(
        "S1_SOP_d0", "S1_SOP_d0", "S2_2Ab_d0", "S3_no2Ab_d0",
        "S4_SOP_d3"
      ), log2_ratio = c(NA, 1, 0, 0, NA)),
      row.names = c(NA, -5L), class = c("tbl_df", "tbl", "data.frame")
    )

  expect_no_warning(
    p <-
      plot_data %>%
      plot_violin(
        x = "sample_alias",
        y = "log2_ratio",
        fill = "sample_alias",
        y_label = "Log2 ratio Proximity Score",
        summarize = FALSE,
        palette = c("red", "blue", "green", "orange"),
        use_jitter = TRUE,
        use_grid = TRUE,
        hline = 0
      )
  )

  jitter_position <- p$layers[[length(p$layers)]]$position
  expect_s3_class(jitter_position, "PositionQuasirandom")
  expect_equal(jitter_position$dodge.width, 0.9)

  # When all displayed_cell_types are present as factor levels on facet_var,
  # empty facets should still be shown without emitting min/max Inf warnings
  sparse_facet_data <-
    tibble(
      sample_alias = factor(rep(c("S1", "S2"), each = 10), levels = c("S1", "S2", "S3")),
      normcount = seq_len(20),
      condition = factor(rep(c("A", "B"), each = 10)),
      celltype = factor(rep(c("T", "B"), each = 10), levels = displayed_cell_types)
    ) %>%
    complete(
      sample_alias = levels(sample_alias),
      celltype = displayed_cell_types
    ) %>%
    mutate(celltype = factor(celltype, displayed_cell_types))

  p_facet <- plot_violin(
    sparse_facet_data,
    x = "sample_alias",
    y = "normcount",
    fill = "condition",
    summarize = FALSE,
    palette = c(A = "red", B = "blue"),
    facet_var = "celltype",
    use_jitter = TRUE,
    use_grid = TRUE
  )

  w_facet <- testthat::capture_warnings(built <- ggplot2::ggplot_build(p_facet))
  expect_false(any(grepl("no non-missing arguments to (min|max)", w_facet)))
  expect_false(any(grepl("stat_ydensity", w_facet)))
  expect_equal(as.character(built$layout$layout$celltype), displayed_cell_types)
  expect_equal(built$layout$panel_params[[1]]$x$get_limits(), c("S1", "S2", "S3"))
})

test_that("PSOCK parallel rendering in tabset_plotlist matches sequential output", {
  skip_on_ci()
  on.exit(pixelatorES:::stop_psock_cluster(), add = TRUE)

  pg_data <- get_test_data(type = "default") |>
    subset(features = rownames(get_test_data(type = "default"))[1:3])
  set.seed(101)
  pg_data[["celltype"]] <- sample(c("B", "T"), ncol(pg_data), replace = TRUE)
  pg_data[["seurat_clusters"]] <- factor(1)
  pg_data[["condition"]] <- "good"
  sample_sheet <- read_samplesheet(test_samplesheet(type = "default"))
  proximity_scores <- suppressWarnings(filter_proximity_scores(
    pg_data,
    list(control_markers = c("mIgG1", "mIgG2a", "mIgG2b"))
  ))

  es_seq <- test_es_data(
    samplesheet = sample_sheet,
    pxl_data_processed = pg_data,
    proximity = proximity_scores,
    params = modifyList(default_params, list(mc_cores = 1L))
  )
  es_par <- test_es_data(
    samplesheet = sample_sheet,
    pxl_data_processed = pg_data,
    proximity = proximity_scores,
    params = modifyList(default_params, list(mc_cores = 2L))
  )
  sample_palette <- create_sample_palette(es_seq$samplesheet$condition)

  dir_seq <- tempfile("knit_seq_")
  dir_par <- tempfile("knit_par_")
  dir.create(dir_seq, recursive = TRUE)
  dir.create(dir_par, recursive = TRUE)
  on.exit(unlink(c(dir_seq, dir_par), recursive = TRUE), add = TRUE)

  env_seq <- new.env(parent = globalenv())
  env_seq$es_data <- es_seq
  env_seq$sample_palette <- sample_palette

  env_par <- new.env(parent = globalenv())
  env_par$es_data <- es_par
  env_par$sample_palette <- sample_palette

  rmd_lines <- c(
    "```{r setup, include=FALSE}",
    "knitr::knit_hooks$set(plot = convert_png_to_webp)",
    "```",
    "```{r abundance_per_marker, fig.height=8, results='asis', echo=FALSE}",
    "plots <- component_abundance_per_marker(es_data, sample_palette, per_celltype = TRUE)",
    "tabset_plotlist(plots, level = 5)",
    "```",
    "```{r proximity_per_marker, fig.height=8, results='asis', echo=FALSE}",
    "plots_prox <- component_proximity_per_marker(es_data, sample_palette, per_celltype = TRUE)",
    "tabset_plotlist(plots_prox, level = 5)",
    "```"
  )
  writeLines(rmd_lines, file.path(dir_seq, "test.Rmd"))
  writeLines(rmd_lines, file.path(dir_par, "test.Rmd"))

  old_wd <- setwd(dir_seq)
  on.exit(setwd(old_wd), add = TRUE)
  knitr::knit("test.Rmd", output = "test.md", quiet = TRUE, envir = env_seq)
  setwd(dir_par)
  knitr::knit("test.Rmd", output = "test.md", quiet = TRUE, envir = env_par)
  setwd(old_wd)

  expect_equal(
    readLines(file.path(dir_seq, "test.md")),
    readLines(file.path(dir_par, "test.md"))
  )

  webp_seq <- sort(list.files(file.path(dir_seq, "figure"), pattern = "\\.webp$", full.names = TRUE))
  webp_par <- sort(list.files(file.path(dir_par, "figure"), pattern = "\\.webp$", full.names = TRUE))
  expect_equal(basename(webp_seq), basename(webp_par))
  expect_length(webp_par, 6L)

  # When knitr.in.progress is TRUE and mc_cores > 1, main-thread plot construction
  # is deferred to the worker pool (elements are lightweight ggplot placeholders)
  expect_true(inherits(env_par$plots[[1]], "ggplot"))
  expect_false(inherits(env_par$plots[[1]], "patchwork"))
  expect_true(inherits(env_seq$plots[[1]], "patchwork"))
  expect_true(inherits(env_par$plots_prox[[1]], "ggplot"))
  expect_false(inherits(env_par$plots_prox[[1]], "patchwork"))
  expect_true(inherits(env_seq$plots_prox[[1]], "patchwork"))

  # Verify sequential fallback in title_plotlist rebuilds deferred plots from spec
  # when parallel rendering does not run (e.g., mc_cores = 1 override)
  p_rebuilt <- pixelatorES:::.build_single_plot_from_spec(
    attr(env_par$plots, "parallel_render_spec"),
    1L
  )
  expect_true(inherits(p_rebuilt, "patchwork"))
})


