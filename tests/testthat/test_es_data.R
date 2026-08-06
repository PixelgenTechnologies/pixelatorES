test_that("es_data construction works as expected", {
  object <- pixelatorES:::new_es_data(list())

  expect_s3_class(object, "es_data")
  expect_equal(
    object[setdiff(names(object), "extractors")],
    list(
      params = list(workflow = "amplicon_demux"),
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

  expect_error(
    pixelatorES:::new_es_data(list(workflow = "unknown"))
  )

  params <- list(sample_sheet = test_samplesheet())
  expect_no_error(built <- build_es_data(params))
  expect_equal(
    built[c("params", "samplesheet", "diagnostics")],
    list(
      params = list(
        sample_sheet = test_samplesheet(),
        workflow = "amplicon_demux"
      ),
      samplesheet = structure(
        list(
          sample = c("S01_PBMC_unstimulated_S1", "S02_PHA_S2"),
          sample_alias = c("S1", "S2"),
          condition = c("PBMC", "PHA")
        ),
        row.names = c(NA, -2L),
        class = c("tbl_df", "tbl", "data.frame")
      ),
      diagnostics = list()
    )
  )
  expect_error(
    build_es_data(list(sample_sheet = tempfile(fileext = ".csv")))
  )
})

test_that("Running es_data extractors works as expected", {
  object <- pixelatorES:::new_es_data(list())
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
  object <- pixelatorES:::new_es_data(list())
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
  register_es_data_workflow(
    "test_workflow",
    function() list(pxl_data = identity),
    overwrite = TRUE
  )

  expect_equal(
    list_es_data_workflows(),
    c("amplicon_demux", "test_workflow")
  )
  expect_equal(
    pixelatorES:::new_es_data(
      list(workflow = "test_workflow")
    )$extractors,
    list(pxl_data = identity)
  )
  expect_error(
    register_es_data_workflow(
      "test_workflow",
      function() list()
    )
  )

  register_es_data_workflow(
    "test_workflow",
    function() list(proximity = identity),
    overwrite = TRUE
  )
  expect_equal(
    pixelatorES:::new_es_data(
      list(workflow = "test_workflow")
    )$extractors,
    list(proximity = identity)
  )
})
