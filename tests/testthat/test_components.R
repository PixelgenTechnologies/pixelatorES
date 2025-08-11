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
    ) %>%
    filter(as.character(marker_1) == as.character(marker_2)) %>%
    group_by(marker_1)
  expect_no_error(
    component <- component_proximity_per_marker(
      proximity_scores,
      proximity_score = "log2_ratio",
      sample_palette = c("red", "black")
    )
  )
  expect_s3_class(component[[1]], "ggplot")

  # component_proximity_selected
  expect_no_error(
    component <- component_proximity_selected(
      proximity_scores,
      sample_palette = c("red", "black"),
      proximity_score = "log2_ratio"
    )
  )

  expect_s3_class(component[[1]], "ggplot")

  # component_dimred_plots
  expect_no_error(
    component <- component_dimred_plots(pg_data, sample_palette = c("red", "black", "blue"))
  )

  expect_s3_class(component$PCA, "ggplot")
})
