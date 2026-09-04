.copy_es_data_test_folder <- function(type) {
  src <- test_data_folder(type = type)
  dest <- tempfile(paste0("pixelatorES_test_", type, "_"))
  dir.create(dest)
  copied <- file.copy(
    list.files(src, full.names = TRUE),
    dest,
    recursive = TRUE
  )
  if (!all(copied)) {
    stop("Failed to copy es_data fixture data", call. = FALSE)
  }

  return(dest)
}

.test_workflow_stages <- function() {
  return(list(
    all = c("graph", "layout"),
    pool = "graph",
    pxl_preference = c("graph", "layout")
  ))
}

test_that("es_data construction works as expected", {
  object <- pixelatorES:::new_es_data(list())

  expect_s3_class(object, "es_data")
  expect_equal(object$stages, get_es_workflow_stages("amplicon_demux"))
  expect_equal(
    object[setdiff(names(object), c("stages", "extractors"))],
    list(
      params = list(workflow = "amplicon_demux"),
      diagnostics = list(),
      samplesheet = NULL,
      sample_aliases = NULL,
      effective_samplesheet = NULL,
      file_paths = NULL,
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
  object <- pixelatorES:::new_es_data(params)
  object$samplesheet <- object$extractors$samplesheet(object)

  expect_equal(
    object[c(
      "params",
      "samplesheet",
      "diagnostics"
    )],
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

test_that("es_data printing works as expected", {
  object <- pixelatorES:::new_es_data(list())
  object$samplesheet <- tibble(
    sample = c("sample_1", "sample_2"),
    sample_alias = c("S1", "S2"),
    condition = c("A", "B")
  )
  object$sample_aliases <- c(sample_1 = "S1", sample_2 = "S2")
  object$effective_samplesheet <- object$samplesheet[1, ]
  object$file_paths <- list(
    data_files = tibble(
      sample_alias = c("S1", "S2"),
      filename = c("sample_1.pxl", "sample_2.pxl")
    ),
    qc_files = tibble(
      sample_alias = c("S1", "S2"),
      filename = c("sample_1.report.json", "sample_2.report.json")
    ),
    pool_qc_files = NULL
  )
  object$pxl_data <- suppressWarnings(
    SeuratObject::CreateSeuratObject(
      counts = matrix(
        1,
        nrow = 2,
        ncol = 3,
        dimnames = list(
          c("marker-1", "marker-2"),
          c("cell-1", "cell-2", "cell-3")
        )
      )
    )
  )
  object$qc_raw <- list(
    qc_files = list(S1 = list(), S2 = list()),
    pool_qc_files = NULL
  )
  object$qc <- list(
    read_stats = tibble(sample_alias = c("S1", "S2")),
    seq_saturation = tibble(sample_alias = c("S1", "S2")),
    denoising = NULL
  )
  object$pxl_data_processed <- object$pxl_data
  object$proximity <- tibble(
    sample_alias = rep(c("S1", "S2"), each = 2),
    score = seq_len(4)
  )
  object <- add_es_data_diagnostic(
    object,
    type = "qc_load",
    target = "S2",
    message = "No QC files were found."
  )
  print_output <- local({
    old_options <- options(cli.width = 80, cli.unicode = FALSE)
    on.exit(options(old_options))
    return(capture.output(print(object), type = "message"))
  })

  expect_equal(
    print_output,
    c(
      "",
      "-- es_data ---------------------------------------------------------------------",
      'Workflow: "amplicon_demux"',
      "* params: 1 parameter",
      "* diagnostics: 1 diagnostic",
      "* samplesheet: tibble with 2 samples and 3 columns",
      "* sample_aliases: named character vector with 2 aliases",
      "* effective_samplesheet: tibble with 1 sample and 3 columns",
      "* file_paths: list with 4 discovered files",
      "* pxl_data: Seurat object with 3 cells and 2 features",
      "* qc_raw: raw QC data for 2 samples and 0 pools",
      "* qc: list with 2 formatted metrics",
      "* pxl_data_processed: Seurat object with 3 cells and 2 features",
      "* proximity: tibble with 4 proximity scores",
      "* stages: workflow vocabulary with 9 stages",
      "* extractors: workflow registry with 15 extractors"
    )
  )
})

test_that("Sample aliases work as expected", {
  expect_equal(
    pixelatorES:::.sample_aliases_from_samplesheet(
      tibble(
        sample = c("sample_1", "sample_2"),
        sample_alias = c("S1", "S2"),
        condition = c("A", "B")
      )
    ),
    c(sample_1 = "S1", sample_2 = "S2")
  )

  expect_equal(
    pixelatorES:::.sample_aliases_from_samplesheet(
      tibble(
        pool = c("pool_1", "pool_1", "pool_2"),
        sample = c("sample_1", "sample_2", "sample_3"),
        sample_alias = c("S1", "S2", "S3"),
        condition = c("A", "B", "A")
      )
    ),
    c(
      sample_1 = "S1",
      sample_2 = "S2",
      sample_3 = "S3",
      pool_1 = "pool_1",
      pool_2 = "pool_2"
    )
  )
})

test_that("test_es_data fixtures work as expected", {
  expect_s3_class(test_es_data(), "es_data")

  prox <- data.frame(marker = "CD3")
  expect_equal(test_es_data(proximity = prox)$proximity, prox)

  qc <- list(read_stats = data.frame(n = 1))
  expect_equal(test_es_data(qc = qc)$qc, qc)

  samplesheet <- tibble(
    sample = c("sample_1", "sample_2"),
    sample_alias = c("S1", "S2"),
    condition = c("A", "B")
  )
  object <- test_es_data(samplesheet = samplesheet)
  expect_equal(
    object$sample_aliases,
    c(sample_1 = "S1", sample_2 = "S2")
  )
  expect_equal(object$effective_samplesheet, samplesheet)

  expect_error(test_es_data(not_a_slot = 1))
})

test_that("Samplesheet extractors work as expected", {
  samplesheet <- tibble(
    sample = c("sample_1", "sample_2"),
    sample_alias = c("S1", "S2"),
    condition = c("A", "B")
  )
  counts <- matrix(
    c(1, 0, 0, 1),
    nrow = 2,
    dimnames = list(
      "marker" = c("CD3", "CD4"),
      "cell" = c("cell_1", "cell_2")
    )
  ) %>%
    as("dgCMatrix")
  pxl_data <- CreateSeuratObject(counts)
  pxl_data$sample_alias <- c("S2", "S2")

  object <- pixelatorES:::new_es_data(list())
  object$samplesheet <- samplesheet
  object$extractors <- list(
    pxl_data = function(object) return(pxl_data),
    effective_samplesheet = pixelatorES:::.extract_effective_samplesheet
  )
  object <- pixelatorES:::run_es_data_extractors(object)

  expect_equal(
    object$effective_samplesheet,
    structure(
      list(
        sample = "sample_2",
        sample_alias = "S2",
        condition = "B"
      ),
      row.names = 1L,
      class = c("tbl_df", "tbl", "data.frame")
    )
  )
})

test_that("Running es_data extractors works as expected", {
  object <- pixelatorES:::new_es_data(list())
  object$extractors <- list(
    pxl_data = function(object) return("raw pxl"),
    qc = list(
      read_stats = function(object) {
        return(data.frame(reads = nchar(object$pxl_data)))
      },
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
  report <- list(
    preamble = c("shared/preprocessing.qmd"),
    sections = list(
      list(id = "samples", title = "Samples", child = "shared/samples.qmd")
    )
  )

  register_es_data_workflow(
    "test_workflow",
    function() return(list(pxl_data = identity)),
    report = function() return(report),
    stages = .test_workflow_stages,
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
  expect_equal(get_es_workflow_report("test_workflow"), report)
  expect_error(
    register_es_data_workflow(
      "test_workflow",
      function() return(list()),
      report = function() return(report),
      stages = .test_workflow_stages
    )
  )

  register_es_data_workflow(
    "test_workflow",
    function() return(list(proximity = identity)),
    report = function() return(report),
    stages = .test_workflow_stages,
    overwrite = TRUE
  )
  expect_equal(
    pixelatorES:::new_es_data(
      list(workflow = "test_workflow")
    )$extractors,
    list(proximity = identity)
  )
})

test_that("Workflow stage vocabularies work as expected", {
  amplicon_vocabulary <- get_es_workflow_stages("amplicon_demux")

  expect_equal(
    amplicon_vocabulary,
    list(
      all = c(
        "amplicon", "collapse",
        "demux", "denoise",
        "graph", "analysis",
        "sample_calling",
        "post_analysis", "layout"
      ),
      pool = c("amplicon", "collapse", "demux", "graph"),
      pxl_preference = c("graph", "analysis", "post_analysis", "layout")
    )
  )
  expect_equal(pixelatorES:::new_es_data(list())$stages, amplicon_vocabulary)
  expect_error(get_es_workflow_stages("unknown_workflow"))

  report <- list(
    preamble = c("shared/preprocessing.qmd"),
    sections = list(
      list(id = "samples", title = "Samples", child = "shared/samples.qmd")
    )
  )

  register_es_data_workflow(
    "test_stages_workflow",
    function() return(list(pxl_data = identity)),
    report = function() return(report),
    stages = .test_workflow_stages,
    overwrite = TRUE
  )
  expect_equal(
    get_es_workflow_stages("test_stages_workflow"),
    .test_workflow_stages()
  )
  expect_equal(
    pixelatorES:::new_es_data(
      list(workflow = "test_stages_workflow")
    )$stages,
    .test_workflow_stages()
  )

  expect_error(
    register_es_data_workflow(
      "missing_stages_workflow",
      function() return(list()),
      report = function() return(report),
      overwrite = TRUE
    )
  )
  expect_error(
    register_es_data_workflow(
      "bad_stages_factory_workflow",
      function() return(list()),
      report = function() return(report),
      stages = .test_workflow_stages(),
      overwrite = TRUE
    )
  )

  register_bad_stages <- function(stages) {
    return(register_es_data_workflow(
      "bad_stages_workflow",
      function() return(list()),
      report = function() return(report),
      stages = function() return(stages),
      overwrite = TRUE
    ))
  }

  expect_error(register_bad_stages(list(all = "graph", pool = "graph")))
  expect_error(
    register_bad_stages(c(.test_workflow_stages(), list(extra = "graph")))
  )
  expect_error(register_bad_stages(unname(.test_workflow_stages())))
  expect_error(
    register_bad_stages(
      modifyList(.test_workflow_stages(), list(all = character()))
    )
  )
  expect_error(
    register_bad_stages(
      modifyList(.test_workflow_stages(), list(all = c("graph", "graph")))
    )
  )
  expect_error(
    register_bad_stages(
      modifyList(.test_workflow_stages(), list(all = c("graph", "")))
    )
  )
  expect_error(
    register_bad_stages(
      modifyList(.test_workflow_stages(), list(pool = "denoise"))
    )
  )
  expect_error(
    register_bad_stages(
      modifyList(.test_workflow_stages(), list(pxl_preference = "denoise"))
    )
  )
  expect_no_error(
    register_bad_stages(
      modifyList(.test_workflow_stages(), list(pool = character()))
    )
  )
})

test_that("File discovery uses the object's stage vocabulary", {
  object <- pixelatorES:::new_es_data(list(
    sample_sheet = test_samplesheet(),
    data_folder = test_data_folder(type = "default")
  ))
  object$samplesheet <- object$extractors$samplesheet(object)

  expect_equal(
    nrow(pixelatorES:::.fill_es_data_file_paths(object)$file_paths$data_files),
    2
  )

  # The fixtures only ship layout PXL files, so narrowing the preference to
  # stages the workflow never emitted must leave no data files behind.
  object$stages$pxl_preference <- c("graph", "analysis")
  narrowed <- pixelatorES:::.fill_es_data_file_paths(object)

  expect_equal(nrow(narrowed$file_paths$data_files), 0)
  expect_gt(nrow(narrowed$file_paths$qc_files), 0)
})

test_that("Workflow report recipes work as expected", {
  report <- list(
    preamble = c("shared/preprocessing.qmd"),
    sections = list(
      list(
        id = "samples",
        title = "Samples",
        child = "shared/samples.qmd"
      ),
      list(
        id = "quality_metrics",
        title = "Quality metrics",
        child = "workflows/amplicon_demux/quality_metrics.qmd"
      )
    )
  )

  register_es_data_workflow(
    "test_report_workflow",
    function() return(list(pxl_data = identity)),
    report = function() return(report),
    stages = .test_workflow_stages,
    overwrite = TRUE
  )
  expect_equal(
    get_es_workflow_report("test_report_workflow"),
    report
  )
  expect_equal(
    get_es_workflow_report("amplicon_demux"),
    list(
      preamble = c("shared/preprocessing.qmd"),
      sections = list(
        list(id = "samples", title = "Samples", child = "shared/samples.qmd"),
        list(
          id = "quality_metrics",
          title = "Quality metrics",
          child = "workflows/amplicon_demux/quality_metrics.qmd"
        ),
        list(
          id = "cell_annotation",
          title = "Cell annotation",
          child = "workflows/amplicon_demux/cell_annotation.qmd"
        ),
        list(
          id = "abundance",
          title = "Abundance",
          child = "workflows/amplicon_demux/abundance.qmd"
        ),
        list(
          id = "spatial",
          title = "Spatial metrics",
          child = "workflows/amplicon_demux/spatial.qmd"
        ),
        list(
          id = "run_info",
          title = "Run info",
          child = "shared/run_info.qmd"
        )
      )
    )
  )

  amplicon_report <- get_es_workflow_report("amplicon_demux")
  amplicon_paths <- c(
    amplicon_report$preamble,
    vapply(amplicon_report$sections, function(section) {
      return(section$child)
    }, character(1))
  )
  quarto_root <- system.file("quarto", package = "pixelatorES")
  expect_true(nzchar(quarto_root))
  expect_equal(
    file.exists(file.path(quarto_root, amplicon_paths)),
    rep(TRUE, length(amplicon_paths))
  )

  expect_error(
    register_es_data_workflow(
      "missing_path_workflow",
      function() return(list()),
      report = function() {
        return(list(
          preamble = "shared/preprocessing.qmd",
          sections = list(
            list(
              id = "samples",
              title = "Samples",
              child = "shared/does_not_exist.qmd"
            )
          )
        ))
      },
      stages = .test_workflow_stages,
      overwrite = TRUE
    )
  )

  expect_error(
    register_es_data_workflow(
      "bad_report_workflow",
      function() return(list()),
      report = function() {
        return(list(sections = list(
          list(id = "samples", title = "Samples", child = "samples.qmd")
        )))
      },
      stages = .test_workflow_stages,
      overwrite = TRUE
    )
  )
  expect_error(
    register_es_data_workflow(
      "bad_report_workflow",
      function() return(list()),
      report = function() {
        return(list(
          preamble = character(),
          sections = list(
            list(id = "samples", title = "Samples", child = "samples.qmd")
          )
        ))
      },
      stages = .test_workflow_stages,
      overwrite = TRUE
    )
  )
  expect_error(
    register_es_data_workflow(
      "bad_report_workflow",
      function() return(list()),
      report = function() {
        return(list(
          preamble = "preprocessing.qmd",
          sections = list()
        ))
      },
      stages = .test_workflow_stages,
      overwrite = TRUE
    )
  )
  expect_error(
    register_es_data_workflow(
      "bad_report_workflow",
      function() return(list()),
      report = function() {
        return(list(
          preamble = "preprocessing.qmd",
          sections = list(
            list(id = "samples", title = "Samples")
          )
        ))
      },
      stages = .test_workflow_stages,
      overwrite = TRUE
    )
  )
  expect_error(
    register_es_data_workflow(
      "bad_report_workflow",
      function() return(list()),
      report = function() {
        return(list(
          preamble = "shared/preprocessing.qmd",
          sections = list(
            list(
              id = "samples",
              title = "Samples",
              child = "shared/samples.qmd"
            ),
            list(
              id = "samples",
              title = "Again",
              child = "shared/samples.qmd"
            )
          )
        ))
      },
      stages = .test_workflow_stages,
      overwrite = TRUE
    )
  )
  expect_error(
    register_es_data_workflow(
      "bad_extractors_workflow",
      function() return("not a list"),
      report = function() return(report),
      stages = .test_workflow_stages,
      overwrite = TRUE
    )
  )
  expect_error(
    register_es_data_workflow(
      "bad_report_factory_workflow",
      function() return(list()),
      report = report,
      stages = .test_workflow_stages,
      overwrite = TRUE
    )
  )
  expect_error(get_es_workflow_report("unknown_workflow"))
})

test_that("Partial input failures work as expected", {
  sample_sheet <- read_samplesheet(test_samplesheet())
  data_folder <- .copy_es_data_test_folder("default")

  pxl_file <- list.files(
    data_folder,
    pattern = paste0("^", sample_sheet$sample[[1]], ".*\\.pxl$"),
    recursive = TRUE,
    full.names = TRUE
  )
  unlink(pxl_file)

  qc_file <- list.files(
    data_folder,
    pattern = paste0("^", sample_sheet$sample[[2]], ".*\\.report\\.json$"),
    recursive = TRUE,
    full.names = TRUE
  )[[1]]
  writeLines("{ invalid json", qc_file)

  object <- pixelatorES:::new_es_data(list(
    sample_sheet = test_samplesheet(),
    data_folder = data_folder
  ))
  object$samplesheet <- object$extractors$samplesheet(object)
  object <- pixelatorES:::.fill_es_data_file_paths(object)
  expect_warning(
    object <- pixelatorES:::run_es_data_extractors(
      object,
      object$extractors["pxl_data"]
    ),
    "Skipping PXL"
  )
  object <- pixelatorES:::run_es_data_extractors(
    object,
    object$extractors["effective_samplesheet"]
  )
  expect_warning(
    object <- pixelatorES:::run_es_data_extractors(
      object,
      object$extractors["qc_raw"]
    ),
    "Skipping QC"
  )

  diagnostics <- lapply(
    object$diagnostics,
    function(diagnostic) {
      return(diagnostic[c("type", "target")])
    }
  )

  expect_equal(
    list(
      effective_samplesheet = object$effective_samplesheet,
      qc_samples = names(object$qc_raw$qc_files),
      qc_pools = object$qc_raw$pool_qc_files,
      diagnostics = diagnostics
    ),
    list(
      effective_samplesheet = structure(
        list(
          sample = "S02_PHA_S2",
          sample_alias = "S2",
          condition = "PHA"
        ),
        row.names = 1L,
        class = c("tbl_df", "tbl", "data.frame")
      ),
      qc_samples = "S1",
      qc_pools = NULL,
      diagnostics = list(
        list(type = "pxl_load", target = "S1"),
        list(type = "qc_load", target = "S2")
      )
    )
  )
})

test_that("Duplicate PXL matches soft-fail as expected", {
  sample_sheet <- read_samplesheet(test_samplesheet())
  data_folder <- .copy_es_data_test_folder("default")
  original_pxl <- list.files(
    data_folder,
    pattern = paste0("^", sample_sheet$sample[[1]], ".*\\.pxl$"),
    recursive = TRUE,
    full.names = TRUE
  )
  extra_dir <- file.path(data_folder, "extra")
  dir.create(extra_dir)
  file.copy(original_pxl, file.path(extra_dir, basename(original_pxl)))

  expect_error(
    get_file_paths(
      data_folder = data_folder,
      sample_sheet = sample_sheet,
      stages = get_es_workflow_stages("amplicon_demux")
    )
  )

  expect_warning(
    object <- build_es_data(list(
      sample_sheet = test_samplesheet(),
      data_folder = data_folder,
      test_mode = TRUE,
      workflow = "amplicon_demux"
    )),
    "Multiple PXL files"
  )

  expect_equal(
    list(
      effective_samplesheet = object$effective_samplesheet$sample_alias,
      duplicate_diagnostics = lapply(
        object$diagnostics[vapply(
          object$diagnostics,
          function(diagnostic) {
            return(
              identical(diagnostic$type, "pxl_load") &&
                identical(diagnostic$target, "S1")
            )
          },
          logical(1)
        )],
        function(diagnostic) {
          return(diagnostic[c("type", "target", "message")])
        }
      )
    ),
    list(
      effective_samplesheet = "S2",
      duplicate_diagnostics = list(
        list(
          type = "pxl_load",
          target = "S1",
          message = "Multiple PXL files matched this sample."
        )
      )
    )
  )
})

test_that("QC-only builds work as expected", {
  data_folder <- .copy_es_data_test_folder("default")
  unlink(list.files(
    data_folder,
    pattern = "\\.pxl$",
    recursive = TRUE,
    full.names = TRUE
  ))

  params <- modifyList(
    default_params,
    list(
      sample_sheet = test_samplesheet(),
      data_folder = data_folder,
      test_mode = TRUE
    )
  )
  object <- suppressWarnings(build_es_data(params))

  diagnostics <- lapply(
    object$diagnostics,
    function(diagnostic) {
      return(diagnostic[c("type", "target")])
    }
  )

  expect_equal(
    list(
      pxl_data = object$pxl_data,
      pxl_data_processed = object$pxl_data_processed,
      proximity = object$proximity,
      effective_samples = object$effective_samplesheet$sample_alias,
      qc_raw_samples = names(object$qc_raw$qc_files),
      populated_qc = names(Filter(Negate(is.null), object$qc)),
      diagnostics = diagnostics
    ),
    list(
      pxl_data = NULL,
      pxl_data_processed = NULL,
      proximity = NULL,
      effective_samples = character(),
      qc_raw_samples = c("S1", "S2"),
      populated_qc = c(
        "crossing_edges",
        "degree_distribution",
        "denoising"
      ),
      diagnostics = list(
        list(type = "pxl_load", target = "S1"),
        list(type = "pxl_load", target = "S2"),
        list(type = "extractor", target = "qc$read_stats"),
        list(type = "extractor", target = "qc$sample_hash_stats"),
        list(type = "extractor", target = "qc$seq_saturation"),
        list(type = "extractor", target = "qc$denoising_detail"),
        list(type = "extractor", target = "qc$coreness"),
        list(type = "extractor", target = "qc$top_markers"),
        list(type = "extractor", target = "pxl_data_processed"),
        list(type = "extractor", target = "proximity")
      )
    )
  )
})

test_that("es_data parity with legacy preprocessing works as expected", {
  .legacy_preprocessing <- function(params) {
    sample_sheet <- read_samplesheet(params$sample_sheet)
    sample_aliases <-
      pixelatorES:::.sample_aliases_from_samplesheet(sample_sheet)

    file_paths <- get_file_paths(
      data_folder = params$data_folder,
      sample_sheet = sample_sheet,
      stages = get_es_workflow_stages("amplicon_demux")
    )

    pg_data <-
      load_pxl_data_list(
        params$data_folder,
        file_paths$data_files,
        sample_sheet
      ) %>%
      merge_data(sample_sheet)

    sample_qc_metrics <- read_qc_files(file_paths, sample_sheet)
    qc_metrics_tables <- get_qc_metrics(
      pg_data,
      sample_qc_metrics,
      sample_sheet
    )

    set.seed(37)
    pg_data_processed <- process_data(pg_data, params)
    proximity_scores <- filter_proximity_scores(
      pg_data_processed,
      params,
      sample_levels = sample_aliases
    )
    component_meta <-
      FetchData(
        pg_data_processed,
        vars = c("sample_alias", "condition", "seurat_clusters", "celltype")
      ) %>%
      as_tibble(rownames = "sample_component")
    proximity_scores <- complete_proximity_scores(
      proximity_scores,
      only_self = FALSE,
      component_meta = component_meta
    )
    proximity_scores <- set_sample_levels(proximity_scores, sample_aliases)

    return(list(
      samplesheet = sample_sheet,
      sample_aliases = sample_aliases,
      pxl_data = pg_data,
      qc_raw = sample_qc_metrics,
      qc = qc_metrics_tables,
      pxl_data_processed = pg_data_processed,
      proximity = proximity_scores
    ))
  }

  .seurat_parity_structure <- function(object) {
    meta <-
      object[[]] %>%
      as_tibble(rownames = "cell_id") %>%
      arrange(cell_id)

    reductions <- Reductions(object)
    if (length(reductions) == 0) {
      reduction_embeddings <- list()
    } else {
      reduction_embeddings <- lapply(
        set_names(reductions),
        function(reduction) {
          return(Embeddings(object, reduction)[meta$cell_id, , drop = FALSE])
        }
      )
    }

    return(list(
      counts = as.matrix(LayerData(object, "counts"))[, meta$cell_id, drop = FALSE],
      meta = meta,
      reductions = reduction_embeddings
    ))
  }

  for (data_type in c("default", "hashing")) {
    params <- modifyList(
      default_params,
      list(
        sample_sheet = test_samplesheet(type = data_type),
        data_folder = .copy_es_data_test_folder(data_type),
        test_mode = TRUE,
        workflow = "amplicon_demux"
      )
    )

    legacy <- .legacy_preprocessing(params)

    params$data_folder <- .copy_es_data_test_folder(data_type)
    set.seed(37)
    built <- build_es_data(params)

    expect_equal(
      list(
        samplesheet = built$samplesheet,
        sample_aliases = built$sample_aliases,
        effective_samplesheet = built$effective_samplesheet,
        diagnostics = built$diagnostics,
        qc_raw = built$qc_raw,
        qc = built$qc,
        proximity = built$proximity,
        pxl_data = built$pxl_data,
        pxl_data_processed = .seurat_parity_structure(built$pxl_data_processed)
      ),
      list(
        samplesheet = legacy$samplesheet,
        sample_aliases = legacy$sample_aliases,
        effective_samplesheet = legacy$samplesheet,
        diagnostics = list(),
        qc_raw = legacy$qc_raw,
        qc = legacy$qc,
        proximity = legacy$proximity,
        pxl_data = NULL,
        pxl_data_processed = .seurat_parity_structure(legacy$pxl_data_processed)
      )
    )
  }
})
