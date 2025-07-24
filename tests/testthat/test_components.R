test_that("`component_cell_recovery` works as expected", {
  sample_qc_metrics <- get_test_qc_metrics()

  expect_no_error(component <- component_cell_recovery(sample_qc_metrics, sample_levels = NULL))
})

test_that("`component_node_degree` works as expected", {
  pg_data <- get_test_data()

  expect_no_error(
    component <- component_node_degree(pg_data, sample_levels = NULL)
  )

  expect_s3_class(component$plot, "ggplot")
  expect_s3_class(component$table, "datatables")
})
