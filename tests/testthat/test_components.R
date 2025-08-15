library(dplyr)
pg_data <- get_test_data()
pg_data_small <- get_test_data(concatenate = FALSE)
sample_qc_metrics <- get_test_qc_metrics()
sample_sheet <- read_samplesheet(test_samplesheet())[1, ]
qc_metrics_tables <-
  get_qc_metrics(pg_data, sample_qc_metrics[1], sample_sheet)

test_that("Components work as expected", {
  # component_control_markers
  expect_no_error(
    component <- component_control_markers(pg_data)
  )

  expect_s3_class(component$p1, "ggplot")
  expect_s3_class(component$p2, "ggplot")
  expect_s3_class(component$tabl, "datatables")

  # component_cell_recovery
  expect_no_error(component <- component_cell_recovery(sample_qc_metrics, sample_levels = NULL))

  # component_node_degree
  expect_no_error(
    component <- component_node_degree(pg_data, sample_levels = NULL)
  )

  expect_s3_class(component$plot, "ggplot")
  expect_s3_class(component$table, "datatables")

  # component_node_edge_count
  expect_no_error(
    component <- component_node_edge_count(sample_qc_metrics, sample_levels = NULL)
  )

  for (plot in component$plots) expect_s3_class(plot, "ggplot")
  expect_s3_class(component$table, "datatables")

  # component_sequencing_saturation
  expect_no_error(
    component <- component_sequencing_saturation(qc_metrics_tables, sample_levels = NULL)
  )

  for (plot in component$plots) expect_s3_class(plot, "ggplot")
  expect_s3_class(component$table, "datatables")

  # component_sequencing_reads_per_cell
  expect_no_error(component <- component_sequencing_reads_per_cell(pg_data))
  expect_s3_class(component$plot, "ggplot")
  expect_s3_class(component$table, "datatables")

  # component_abundance_per_celltype
  temp <-
    pg_data %>%
    subset(features = rownames(pg_data)[1:3])

  temp[["l1_annotation_summary"]] <- sample(
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

  # component_proximity_per_marker
  temp <-
    pg_data_small %>%
    subset(features = rownames(pg_data_small)[1:3])

  temp[["l1_annotation_summary"]] <- sample(
    c("CD4 T", "CD8 T", "B"),
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
  expect_equal(
    component[[1]]$data,
    structure(list(
      sample_alias = c(
        "S1", "S1", "S1", "S1", "S1",
        "S1"
      ), l1_annotation_summary = c(
        "B", "CD4 T", "CD4 T", "CD8 T",
        "Mono", "NK"
      ), marker_1 = structure(c(1L, 1L, 1L, NA, NA, NA), levels = c(
        "CD11b",
        "B2M", "HLA-ABC"
      ), class = "factor"), marker_2 = structure(c(
        1L,
        1L, 1L, NA, NA, NA
      ), levels = c("CD11b", "B2M", "HLA-ABC"), class = "factor"),
      join_count = c(0, 0, 57, NA, NA, NA), join_count_expected_mean = c(
        0,
        0.03, 45.18, NA, NA, NA
      ), join_count_expected_sd = c(
        0, 0.171446607997765,
        6.58338895258438, NA, NA, NA
      ), join_count_z = c(
        0, -0.03,
        1.79542786931341, NA, NA, NA
      ), join_count_p = c(
        0.5, 0.488033526585887,
        0.0362927777357692, NA, NA, NA
      ), log2_ratio = c(
        0, 0, 0.335277648546382,
        NA, NA, NA
      ), sample_component = c(
        "2708240b908e2eba", "c3c393e9a17c1981",
        "0a45497c6bfbfb22", NA, NA, NA
      ), count_1 = c(
        12L, 17L, 929L,
        NA, NA, NA
      ), count_2 = c(12L, 17L, 929L, NA, NA, NA), p1 = c(
        0.00216723857684667,
        0.0021783700666325, 0.312163978494624, NA, NA, NA
      ), p2 = c(
        0.00216723857684667,
        0.0021783700666325, 0.312163978494624, NA, NA, NA
      ), condition = c(
        "good",
        "good", "good", NA, NA, NA
      ), seurat_clusters = c(
        "1", "1",
        "1", NA, NA, NA
      )
    ), class = c("tbl_df", "tbl", "data.frame"), row.names = c(NA, -6L))
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
  expect_equal(
    component[[1]]$data,
    structure(list(
      sample_alias = c(
        "S1", "S1", "S1", "S1", "S1",
        "S1"
      ), l1_annotation_summary = c(
        "B", "CD4 T", "CD4 T", "CD8 T",
        "Mono", "NK"
      ), marker_1 = structure(c(1L, 1L, 1L, NA, NA, NA), levels = c(
        "CD11b",
        "B2M", "HLA-ABC"
      ), class = "factor"), marker_2 = structure(c(
        1L,
        1L, 1L, NA, NA, NA
      ), levels = c("CD11b", "B2M", "HLA-ABC"), class = "factor"),
      join_count = c(0, 0, 57, NA, NA, NA), join_count_expected_mean = c(
        0,
        0.03, 45.18, NA, NA, NA
      ), join_count_expected_sd = c(
        0, 0.171446607997765,
        6.58338895258438, NA, NA, NA
      ), join_count_z = c(
        0, -0.03,
        1.79542786931341, NA, NA, NA
      ), join_count_p = c(
        0.5, 0.488033526585887,
        0.0362927777357692, NA, NA, NA
      ), log2_ratio = c(
        0, 0, 0.335277648546382,
        NA, NA, NA
      ), sample_component = c(
        "2708240b908e2eba", "c3c393e9a17c1981",
        "0a45497c6bfbfb22", NA, NA, NA
      ), count_1 = c(
        12L, 17L, 929L,
        NA, NA, NA
      ), count_2 = c(12L, 17L, 929L, NA, NA, NA), p1 = c(
        0.00216723857684667,
        0.0021783700666325, 0.312163978494624, NA, NA, NA
      ), p2 = c(
        0.00216723857684667,
        0.0021783700666325, 0.312163978494624, NA, NA, NA
      ), condition = c(
        "good",
        "good", "good", NA, NA, NA
      ), seurat_clusters = c(
        "1", "1",
        "1", NA, NA, NA
      )
    ), class = c("tbl_df", "tbl", "data.frame"), row.names = c(NA, -6L))
  )

  # component_dimred_plots
  expect_no_error(
    component <- component_dimred_plots(pg_data, sample_palette = c("red", "black", "blue"))
  )

  expect_s3_class(component$PCA, "ggplot")
})
