seur <-
  minimal_pna_pxl_file() %>%
  ReadPNA_Seurat(load_proximity_scores = FALSE)

seur <-
  merge(
    seur,
    y = list(seur, seur, seur, seur, seur)
  ) %>%
  JoinLayers()

seur[[]]$sample_alias <-
  "S1"

seur <-
  seur %>%
  NormalizeData() %>%
  ScaleData() %>%
  RunPCA(npcs = 2)

test_that("Doublet Detection component works as expected", {
  expect_no_error(
    component <- component_doublet_detection(
      seur,
      c("red", "black")
    )
  )

  expect_s3_class(component$barplot, "ggplot")
  expect_s3_class(component$dimred_plot, "ggplot")
})
