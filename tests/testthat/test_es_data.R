test_that("new_es_data creates the stable data shape with workflow extractors", {
  samplesheet <- tibble::tibble(
    sample = "S1",
    sample_alias = "S1",
    condition = "A"
  )
  object <- pixelatorES:::new_es_data(
    samplesheet = samplesheet,
    meta = list(run_id = "test-run")
  )

  expect_s3_class(object, "es_data")
  expect_named(
    object,
    c(
      "meta",
      "diagnostics",
      "samplesheet",
      "effective_samplesheet",
      "pxl_data",
      "qc_raw",
      "qc",
      "pxl_data_processed",
      "proximity",
      "extractors"
    )
  )
  expect_identical(object$meta$workflow, "amplicon_demux")
  expect_identical(object$meta$run_id, "test-run")
  expect_identical(object$samplesheet, samplesheet)
  expect_null(object$effective_samplesheet)
  expect_length(object$diagnostics, 0)
  expect_null(object$pxl_data)
  expect_length(object$qc, 0)
  expect_named(
    object$extractors,
    c("pxl_data", "qc_raw", "qc", "pxl_data_processed", "proximity")
  )
  expect_true(is.function(object$extractors$pxl_data))
  expect_true(is.list(object$extractors$qc))
  expect_true(is.function(object$extractors$qc$crossing_edges))
})

test_that("new_es_data rejects unknown workflows", {
  expect_error(pixelatorES:::new_es_data(workflow = "unknown"))
})

test_that("build_es_data runs the default workflow", {
  object <- build_es_data()

  expect_s3_class(object, "es_data")
  expect_identical(object$meta$workflow, "amplicon_demux")
  expect_null(object$samplesheet)
  expect_null(object$effective_samplesheet)
  expect_null(object$pxl_data)
  expect_null(object$qc$read_stats)
  expect_length(object$diagnostics, 0)
})

test_that("extractors fill nested slots in registry order", {
  object <- pixelatorES:::new_es_data()
  object$extractors <- list(
    pxl_data = function(object) "raw pxl",
    qc = list(
      read_stats = function(object) {
        expect_identical(object$pxl_data, "raw pxl")
        data.frame(sample_alias = "sample-1", reads = 10)
      }
    ),
    pxl_data_processed = function(object) paste(object$pxl_data, "processed")
  )

  object <- pixelatorES:::run_es_data_extractors(object)

  expect_identical(object$pxl_data, "raw pxl")
  expect_identical(object$pxl_data_processed, "raw pxl processed")
  expect_identical(object$qc$read_stats$reads, 10)
  expect_length(object$diagnostics, 0)
})

test_that("extractor failures are isolated and recorded", {
  object <- pixelatorES:::new_es_data()
  object$extractors <- list(
    proximity = function(object) stop("proximity unavailable"),
    qc = list(
      read_stats = function(object) data.frame(reads = 10)
    )
  )

  object <- pixelatorES:::run_es_data_extractors(object)

  expect_null(object$proximity)
  expect_identical(object$qc$read_stats$reads, 10)
  expect_length(object$diagnostics, 1)
  expect_identical(object$diagnostics[[1]]$type, "extractor")
  expect_identical(object$diagnostics[[1]]$target, "proximity")
  expect_match(object$diagnostics[[1]]$message, "proximity unavailable")
})

test_that("nested extractor failures use a dotted slot target", {
  object <- pixelatorES:::new_es_data()
  object$extractors <- list(
    qc = list(
      crossing_edges = function(object) stop("missing edges")
    )
  )

  object <- pixelatorES:::run_es_data_extractors(object)

  expect_null(object$qc$crossing_edges)
  expect_identical(object$diagnostics[[1]]$target, "qc$crossing_edges")
})

test_that("add_es_data_diagnostic records type target and message", {
  object <- pixelatorES:::new_es_data()
  object <- pixelatorES:::add_es_data_diagnostic(
    object,
    type = "pxl_load",
    target = "S1",
    message = "corrupt pxl"
  )

  expect_length(object$diagnostics, 1)
  expect_identical(
    object$diagnostics[[1]],
    list(type = "pxl_load", target = "S1", message = "corrupt pxl")
  )
  expect_error(
    pixelatorES:::add_es_data_diagnostic(
      object,
      type = "not_a_type",
      target = "S1",
      message = "bad"
    )
  )
})
