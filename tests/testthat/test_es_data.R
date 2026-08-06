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
  register_es_data_workflow(
    "test_workflow",
    function() return(list(pxl_data = identity)),
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
      function() return(list())
    )
  )

  register_es_data_workflow(
    "test_workflow",
    function() return(list(proximity = identity)),
    overwrite = TRUE
  )
  expect_equal(
    pixelatorES:::new_es_data(
      list(workflow = "test_workflow")
    )$extractors,
    list(proximity = identity)
  )
})

test_that("es_data parity with legacy preprocessing works as expected", {
  .copy_parity_data_folder <- function(type) {
    src <- test_data_folder(type = type)
    dest <- tempfile(paste0("pixelatorES_parity_", type, "_"))
    dir.create(dest)
    copied <- file.copy(
      list.files(src, full.names = TRUE),
      dest,
      recursive = TRUE
    )
    if (!all(copied)) {
      stop("Failed to copy parity fixture data", call. = FALSE)
    }
    return(dest)
  }

  .legacy_preprocessing <- function(params) {
    sample_sheet <- read_samplesheet(params$sample_sheet)
    sample_aliases <-
      sample_sheet %>%
      select(sample, sample_alias) %>%
      deframe()

    if ("pool" %in% names(sample_sheet)) {
      pool_aliases <- unique(sample_sheet$pool)
      sample_aliases <-
        c(sample_aliases, set_names(setdiff(pool_aliases, sample_aliases)))
    }

    file_paths <- get_file_paths(
      data_folder = params$data_folder,
      sample_sheet = sample_sheet
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

    return(list(
      samplesheet = sample_sheet,
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
        data_folder = .copy_parity_data_folder(data_type),
        test_mode = TRUE,
        workflow = "amplicon_demux"
      )
    )

    legacy <- .legacy_preprocessing(params)

    params$data_folder <- .copy_parity_data_folder(data_type)
    set.seed(37)
    built <- build_es_data(params)

    expect_equal(
      list(
        samplesheet = built$samplesheet,
        effective_samplesheet = built$effective_samplesheet,
        diagnostics = built$diagnostics,
        qc_raw = built$qc_raw,
        qc = built$qc,
        proximity = built$proximity,
        pxl_data = .seurat_parity_structure(built$pxl_data),
        pxl_data_processed = .seurat_parity_structure(built$pxl_data_processed)
      ),
      list(
        samplesheet = legacy$samplesheet,
        effective_samplesheet = legacy$samplesheet,
        diagnostics = list(),
        qc_raw = legacy$qc_raw,
        qc = legacy$qc,
        proximity = legacy$proximity,
        pxl_data = .seurat_parity_structure(legacy$pxl_data),
        pxl_data_processed = .seurat_parity_structure(legacy$pxl_data_processed)
      )
    )
  }
})
