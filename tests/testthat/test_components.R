library(dplyr)
pg_data <- get_test_data()
pg_data_small <- get_test_data(concatenate = FALSE)
sample_qc_metrics <- get_test_qc_metrics()
sample_sheet <- read_samplesheet(test_samplesheet())
qc_metrics_tables <-
  get_qc_metrics(pg_data, sample_qc_metrics, sample_sheet)

test_that("Components work as expected", {
  # component_control_markers
  expect_no_error(
    component <- component_control_markers(pg_data)
  )

  expect_s3_class(component$p1, "ggplot")
  expect_no_error(
    ggplot2::ggplot_build(component$p1)
  )
  expect_s3_class(component$p2, "ggplot")
  expect_no_error(
    ggplot2::ggplot_build(component$p2)
  )
  expect_s3_class(component$tabl, "datatables")

  # component_cell_recovery
  expect_no_error(component <- component_cell_recovery(sample_qc_metrics, sample_levels = NULL))
  expect_no_error(
    ggplot2::ggplot_build(component$plot[[1]])
  )

  # component_node_degree
  expect_no_error(
    component <- component_node_degree(pg_data, sample_levels = NULL)
  )

  expect_s3_class(component$plot, "ggplot")
  expect_no_error(
    ggplot2::ggplot_build(component$plot)
  )
  expect_s3_class(component$table, "datatables")

  # component_node_edge_count
  expect_no_error(
    component <- component_node_edge_count(sample_qc_metrics, sample_levels = NULL)
  )

  for (plot in component$plots) {
    expect_s3_class(plot, "ggplot")
    expect_no_error(
      ggplot2::ggplot_build(plot)
    )
  }
  expect_s3_class(component$table, "datatables")

  # component_sequencing_saturation
  expect_no_error(
    component <- component_sequencing_saturation(qc_metrics_tables, sample_levels = NULL)
  )

  for (plot in component$plots) {
    expect_s3_class(plot, "ggplot")
    expect_no_error(
      ggplot2::ggplot_build(plot)
    )
  }
  expect_s3_class(component$table, "datatables")

  # component_sequencing_reads_per_cell
  expect_no_error(component <- component_sequencing_reads_per_cell(pg_data))
  expect_s3_class(component$plot, "ggplot")
  expect_no_error(
    ggplot2::ggplot_build(component$plot)
  )
  expect_s3_class(component$table, "datatables")

  # component_abundance_per_celltype
  temp <-
    pg_data %>%
    subset(features = rownames(pg_data)[1:3])

  temp[["celltype"]] <- sample(
    c("CD4 T", "CD8 T", "B"),
    ncol(temp),
    replace = TRUE
  )
  temp[["seurat_clusters"]] <- 1

  expect_no_error(
    component <- component_abundance_per_celltype(
      temp,
      params = list(
        control_markers = c("mIgG1", "mIgG2a", "mIgG2b")
      ),
      sample_palette = c("red", "black")
    )
  )

  expect_s3_class(component[[1]], "ggplot")
  expect_no_error(
    ggplot2::ggplot_build(component[[1]])
  )
  # component_proximity_per_marker
  temp <-
    pg_data_small %>%
    subset(features = rownames(pg_data_small)[1:3])

  set.seed(37)
  temp[["celltype"]] <- sample(
    c("T", "B"),
    ncol(temp),
    replace = TRUE
  )
  temp[["seurat_clusters"]] <- 1
  temp[["condition"]] <- "good"

  proximity_scores <-
    filter_proximity_scores(
      temp,
      list(
        control_markers = c("mIgG1", "mIgG2a", "mIgG2b")
      )
    )

  expect_no_error(
    component <- component_proximity_per_marker(
      proximity_scores,
      sample_palette = c("red", "black")
    )
  )
  expect_s3_class(component[[1]], "ggplot")
  expect_no_error(
    ggplot2::ggplot_build(component[[1]])
  )
  expect_equal(
    component[[1]]$data,
    structure(list(
      sample_alias = structure(c(1L, 1L, 1L, 1L, 1L), levels = "S1", class = "factor"),
      celltype = c("B", "Mono & DC", "NK", "Platelets", "T"), marker_1 = structure(c(
        1L,
        NA, NA, NA, NA
      ), levels = c("CD11b", "B2M", "HLA-ABC"), class = "factor"),
      marker_2 = structure(c(1L, NA, NA, NA, NA), levels = c(
        "CD11b",
        "B2M", "HLA-ABC"
      ), class = "factor"), join_count = c(
        57,
        NA, NA, NA, NA
      ), join_count_expected_mean = c(
        45.18, NA,
        NA, NA, NA
      ), join_count_expected_sd = c(
        6.58338895258438,
        NA, NA, NA, NA
      ), join_count_z = c(
        1.79542786931341, NA, NA,
        NA, NA
      ), join_count_p = c(
        0.0362927777357692, NA, NA, NA,
        NA
      ), log2_ratio = c(0.335277648546382, NA, NA, NA, NA), sample_component = c(
        "0a45497c6bfbfb22",
        NA, NA, NA, NA
      ), count_1 = c(929L, NA, NA, NA, NA), count_2 = c(
        929L,
        NA, NA, NA, NA
      ), p1 = c(0.312163978494624, NA, NA, NA, NA), p2 = c(0.312163978494624, NA, NA, NA, NA), condition = c(
        "good",
        NA, NA, NA, NA
      ), seurat_clusters = c("1", NA, NA, NA, NA)
    ), row.names = c(
      NA,
      -5L
    ), class = c("tbl_df", "tbl", "data.frame"))
  )

  # component_proximity_selected
  set.seed(37)
  expect_no_error(
    component <- component_proximity_selected(
      temp,
      proximity_scores %>%
        filter(as.character(marker_1) == as.character(marker_2)),
      sample_palette = c("red", "black"),
      selected_contrasts = FALSE,
      proximity_score = "log2_ratio"
    )
  )

  expect_named(component, expected = c("B2M/B2M", "CD11b/CD11b", "HLA-ABC/HLA-ABC"))
  expect_s3_class(component[[1]], "ggplot")
  expect_no_error(
    ggplot2::ggplot_build(component[[1]])
  )


  # component_dimred_plots
  expect_no_error(
    component <- component_dimred_plots(pg_data, sample_palette = c("red", "black", "blue"))
  )

  expect_s3_class(component$PCA, "ggplot")

  # component_proximity_heatmap_sample
  expect_no_error(
    component <- component_proximity_heatmap_sample(
      proximity_scores,
      temp,
      heatmap_gradient = c("red", "cyan"),
      n_markers = 2,
      plot_markers = NULL,
      test_mode = TRUE
    )
  )
  # component_proximity_heatmap_celltype
  expect_no_error(
    component <- component_proximity_heatmap_celltype(
      proximity_scores,
      temp,
      heatmap_gradient = c("red", "cyan"),
      n_markers = 2,
      plot_markers = NULL,
      test_mode = TRUE
    )
  )
})
