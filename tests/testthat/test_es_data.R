test_that("es_data construction works as expected", {
  object <- pixelatorES:::new_es_data()

  expect_s3_class(object, "es_data")
  expect_equal(
    object[setdiff(names(object), "extractors")],
    list(
      meta = list(workflow = "amplicon_demux"),
      diagnostics = list(),
      samplesheet = NULL,
      effective_samplesheet = NULL,
      pxl_data = NULL,
      qc_raw = NULL,
      qc = list(),
      pxl_data_processed = NULL,
      proximity = NULL
    )
  )

  samplesheet <- tibble::tibble(
    sample = "S1",
    sample_alias = "S1",
    condition = "A"
  )
  expect_equal(
    pixelatorES:::new_es_data(samplesheet = samplesheet)$samplesheet,
    samplesheet
  )

  expect_error(pixelatorES:::new_es_data(workflow = "unknown"))

  expect_no_error(built <- build_es_data())
  expect_s3_class(built, "es_data")
  expect_equal(built$diagnostics, list())
})

test_that("Running es_data extractors works as expected", {
  object <- pixelatorES:::new_es_data()
  object$extractors <- list(
    pxl_data = function(object) "raw pxl",
    qc = list(
      read_stats = function(object) data.frame(reads = nchar(object$pxl_data)),
      crossing_edges = function(object) stop("missing edges")
    ),
    proximity = function(object) stop("proximity unavailable")
  )

  object <- pixelatorES:::run_es_data_extractors(object)

  expect_equal(
    object[c("pxl_data", "qc", "proximity", "diagnostics")],
    list(
      pxl_data = "raw pxl",
      qc = list(read_stats = data.frame(reads = 7L), crossing_edges = NULL),
      proximity = NULL,
      diagnostics = list(
        list(
          type = "extractor",
          target = "qc$crossing_edges",
          message = "missing edges"
        ),
        list(
          type = "extractor",
          target = "proximity",
          message = "proximity unavailable"
        )
      )
    )
  )
})

test_that("es_data diagnostics work as expected", {
  object <- pixelatorES:::new_es_data()
  object <- pixelatorES:::add_es_data_diagnostic(
    object,
    type = "pxl_load",
    target = "S1",
    message = "corrupt pxl"
  )

  expect_equal(
    object$diagnostics,
    list(list(type = "pxl_load", target = "S1", message = "corrupt pxl"))
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

test_that("es_data workflow registration works as expected", {
  expect_equal(list_es_data_workflows(), "amplicon_demux")

  register_es_data_workflow(
    "test_workflow",
    function() list(pxl_data = identity)
  )

  expect_equal(
    list_es_data_workflows(),
    c("amplicon_demux", "test_workflow")
  )
  expect_equal(
    pixelatorES:::new_es_data(workflow = "test_workflow")$extractors,
    list(pxl_data = identity)
  )

  # Can't overwrite registered workflow
  expect_error(
    register_es_data_workflow(
      "test_workflow",
      function() list()
    )
  )
})
