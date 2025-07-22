
test_that("`component_cell_recovery` works as expected", {

  sample_qc_metrics <- get_test_qc_metrics()

  expect_no_error(component <- component_cell_recovery(sample_qc_metrics, sample_levels = NULL))

})
