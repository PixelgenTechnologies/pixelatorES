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
    c("::: {.panel-tabset .nav-pills}", "## a", "", "", "", "## b",
      "", "", "", "## c", "", "", "", ":::")
  )
  expect_equal(
    capture.output(tabset_plotlist(plot_list, level = 3)),
    c("::: {.panel-tabset .nav-pills}", "### a", "", "", "", "### b",
      "", "", "", "### c", "", "", "", ":::")
  )


  nested_plot_list <-
    list(
      "a" = g, "b" = g, "c" = g,
      plot_list = plot_list
    )


  expect_equal(
    capture.output(tabset_nested_plotlist(nested_plot_list, level = 3)),
    c("", "", "::: {.panel-tabset .nav-pills}", "### a", "", "",
      "", "### b", "", "", "", "### c", "", "", "", "### plot_list",
      "", "", "::: {.panel-tabset .nav-pills}", "#### a", "", "", "",
      "#### b", "", "", "", "#### c", "", "", "", ":::", ":::")
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
