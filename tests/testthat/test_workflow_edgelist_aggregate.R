test_that("edgelist_aggregate workflow registration and stage vocabulary work as expected", {
  expect_true("edgelist_aggregate" %in% list_es_data_workflows())

  stages <- get_es_workflow_stages("edgelist_aggregate")
  expect_equal(
    stages,
    list(
      all = c(
        "edgelist", "aggregate", "graph", "graph_collapse", "denoise",
        "sample_calling", "analysis", "post_analysis", "layout"
      ),
      pool = c("edgelist", "aggregate", "graph", "graph_collapse", "denoise"),
      pxl_preference = c("graph", "analysis", "post_analysis", "layout")
    )
  )

  recipe <- get_es_workflow_report("edgelist_aggregate")
  expect_equal(
    recipe,
    list(
      preamble = "shared/preprocessing.qmd",
      sections = list(
        list(id = "samples", title = "Samples", child = "shared/samples.qmd"),
        list(
          id = "quality_metrics",
          title = "Quality metrics",
          child = "workflows/edgelist_aggregate/quality_metrics.qmd"
        ),
        list(
          id = "cell_annotation",
          title = "Cell annotation",
          child = "workflows/shared/cell_annotation.qmd"
        ),
        list(
          id = "abundance",
          title = "Abundance",
          child = "workflows/shared/abundance.qmd"
        ),
        list(
          id = "spatial",
          title = "Spatial metrics",
          child = "workflows/shared/spatial.qmd"
        ),
        list(
          id = "run_info",
          title = "Run info",
          child = "shared/run_info.qmd"
        )
      )
    )
  )

  own <- pixelatorES:::.get_es_data_extractors("edgelist_aggregate")
  borrowed <- pixelatorES:::.get_es_data_extractors("amplicon_demux")
  expect_equal(names(own), names(borrowed))
  expect_equal(names(own$qc), c(names(borrowed$qc), "graph_collapse"))
  expect_identical(own$qc_raw, borrowed$qc_raw)
})

test_that("edgelist_aggregate file discovery and multi-FASTQ qc_raw aggregation work as expected", {
  sample_sheet <- read_samplesheet(test_samplesheet("edgelist_aggregate"))
  expect_equal(
    sample_sheet,
    structure(
      list(
        sample = c("sample1", "sample2", "sample3", "sample4"),
        sample_alias = c("S1", "S2", "S3", "S4"),
        condition = c("unstim", "stim", "unstim", "stim"),
        pool = c("pool1", "pool1", "pool2", "pool2")
      ),
      row.names = c(NA, -4L),
      class = c("tbl_df", "tbl", "data.frame")
    )
  )

  run_dir <- materialise_edgelist_aggregate_run()
  on.exit(unlink(run_dir, recursive = TRUE), add = TRUE)

  file_paths <- edgelist_aggregate_file_paths(run_dir)
  expect_equal(sort(file_paths$data_files$sample_alias), paste0("S", 1:4))
  expect_true(all(grepl("\\.layout\\.pxl$", file_paths$data_files$filename)))

  expect_equal(
    unname(file_paths$pool_qc_files$stage[file_paths$pool_qc_files$pool_alias == "pool1"]),
    c("aggregate", "denoise", "edgelist", "edgelist", "graph", "graph_collapse", "sample_calling")
  )
  expect_equal(
    unname(file_paths$pool_qc_files$stage[file_paths$pool_qc_files$pool_alias == "pool2"]),
    c("aggregate", "denoise", "edgelist", "graph", "graph_collapse", "sample_calling")
  )

  qc_raw <- edgelist_aggregate_qc_raw(run_dir)
  expect_equal(
    qc_raw$pool_qc_files$pool1$edgelist,
    list(input_reads = 309597, output_reads = 309597)
  )
  expect_equal(
    qc_raw$pool_qc_files$pool1$edgelist$input_reads,
    qc_raw$pool_qc_files$pool1$aggregate$input_reads
  )
  expect_equal(
    qc_raw$pool_qc_files$pool2$edgelist$input_reads,
    155203
  )
})

test_that("edgelist_aggregate QC extractors, key table, and reads per stage work as expected", {
  run_dir <- materialise_edgelist_aggregate_run()
  on.exit(unlink(run_dir, recursive = TRUE), add = TRUE)

  samplesheet <- read_samplesheet(test_samplesheet("edgelist_aggregate"))
  qc_raw <- edgelist_aggregate_qc_raw(run_dir)
  pxl_data <- get_test_data("edgelist_aggregate")

  object <- test_es_data(
    samplesheet = samplesheet,
    qc_raw = qc_raw,
    pxl_data = pxl_data,
    params = list(workflow = "edgelist_aggregate")
  )

  seq_sat <- pixelatorES:::.extract_edgelist_aggregate_seq_saturation(object)
  expect_equal(
    seq_sat,
    structure(
      list(
        pool = structure(1:2, levels = c("pool1", "pool2"), class = "factor"),
        total_reads = c(309597, 155203),
        deduped_valid_reads = c(44400, 22200),
        valid_reads = c(309597, 155203),
        graph_edges = c(66575L, 84144L),
        graph_proteins = c(40627L, 48052L),
        graph_reads = c(154863L, 203615L),
        valid_reads_saturation = c(85.658775763331, 85.6961527805519),
        graph_edge_saturation = c(57.010389828429, 58.6749502738011),
        graph_node_saturation = c(73.7658446497872, 76.400559880166),
        fraction_valid_reads = c(100, 100),
        fraction_graph_reads = c(50.0208335352087, 131.192696017474)
      ),
      row.names = c(NA, -2L),
      class = c("tbl_df", "tbl", "data.frame")
    )
  )

  deg_dist <- pixelatorES:::.extract_edgelist_aggregate_degree_distribution(object)
  expect_equal(
    head(deg_dist, 4),
    structure(
      list(
        pool = structure(c(1L, 1L, 1L, 1L), levels = c("pool1", "pool2"), class = "factor"),
        umi_type = c("umi1", "umi1", "umi1", "umi1"),
        degree = c(1L, 4L, 5L, 6L),
        n = c(2000L, 7868L, 1882L, 233L),
        percent_nodes = c(8.33333333333333, 32.7833333333333, 7.84166666666667, 0.970833333333333)
      ),
      row.names = c(NA, -4L),
      class = c("tbl_df", "tbl", "data.frame")
    )
  )

  graph_col <- pixelatorES:::.extract_graph_collapse(object)
  expect_equal(
    graph_col,
    structure(
      list(
        pool = structure(1:2, levels = c("pool1", "pool2"), class = "factor"),
        nodes_collapsed_fraction = c(0.435714285714286, 0.129807692307692)
      ),
      row.names = c(NA, -2L),
      class = c("tbl_df", "tbl", "data.frame")
    )
  )

  reads_stage <- component_reads_per_stage(object)
  expect_equal(names(reads_stage$plots), c("pool1", "pool2"))
  expect_equal(
    levels(reads_stage$plots$pool1$data$stage),
    c("fastQ", "edgelist", "graph", "graph_collapse", "denoise", "sample_calling")
  )
})
