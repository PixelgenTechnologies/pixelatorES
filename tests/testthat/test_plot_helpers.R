test_that("Tab setting and title setting work as expected", {
  g <-
    ggplot() +
    theme_void()

  plot_list <-
    list("a" = g, "b" = g, "c" = g)

  expect_equal(
    capture.output(title_plotlist(plot_list)),
    c("## a", "", "", "", "## b", "", "", "", "## c", "", "", "")
  )
  expect_equal(
    capture.output(title_plotlist(plot_list, level = 3)),
    c("### a", "", "", "", "### b", "", "", "", "### c", "", "", "")
  )


  expect_equal(
    capture.output(tabset_plotlist(plot_list)),
    c(
      "::: {.panel-tabset .nav-pills}", "## a", "", "", "", "## b",
      "", "", "", "## c", "", "", "", ":::"
    )
  )
  expect_equal(
    capture.output(tabset_plotlist(plot_list, level = 3)),
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

  expect_no_error(title_plotlist(plot_list))

  ## Throw error in dev mode
  options(pixelatorES.dev_mode = TRUE)

  expect_error(title_plotlist(plot_list))

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
    capture.output(tabset_nested_plotlist(nested_plot_list, level = 3)),
    c(
      "", "", "::: {.panel-tabset .nav-pills}", "### a", "", "",
      "", "### b", "", "", "", "### c", "", "", "", "### plot_list",
      "", "", "::: {.panel-tabset .nav-pills}", "#### a", "", "", "",
      "#### b", "", "", "", "#### c", "", "", "", ":::", ":::"
    )
  )
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
        jitter_alpha = 1,
        hline = 0
      )
  )
})
