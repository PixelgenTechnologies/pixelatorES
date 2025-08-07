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

})
