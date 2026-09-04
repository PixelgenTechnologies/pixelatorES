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

  expect_no_error(
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
})
