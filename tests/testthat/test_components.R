pg_data <- get_test_data()
sample_qc_metrics <- get_test_qc_metrics()

test_that("Components work as expected", {

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
    ncol(temp), replace = TRUE)
  temp[["seurat_clusters"]] <- 1

  expect_no_error(
    component <- component_abundance_per_celltype(
      temp,
      params = list(
        control_markers = c("mIgG1", "mIgG2a", "mIgG2b")
      ),
      sample_palette = c("red", "black"))
  )

  expect_s3_class(component[[1]], "ggplot")

})
