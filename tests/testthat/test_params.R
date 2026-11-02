test_that("Parameter tables work as expected", {
  params <- structure(list(
    sample_sheet = "../../../ES_test_data/test_samplesheet.csv",
    data_folder = "../../../ES_test_data",
    control_markers = "c('mIgG1', 'mIgG2a', 'mIgG2b')",
    molecule_rank_cutoff = 10000L,
    gini_dispersion_cutoff = 0.7,
    control_score_cutoff = 1.5,
    norm_method = "CLR",
    clustering_resolution = 1L,
    run_proximity_analysis = TRUE,
    proximity_count_cutoff = 25L,
    run_qc_metrics = TRUE
  ), class = "knit_param_list")
  sample_sheet <- structure(list(
    sample = c("S01_PBMC_unstimulated_S1", "S02_PHA_S2"),
    design = c("pna-2", "pna-2"),
    panel = c("pna-rnd-158plex-final", "pna-rnd-158plex-final"),
    fastq_1 = c(
      "s3://pixelgen-technologies-ngs/NextSeq2000/250314_VH00725_275_AACNHV5HV/Analysis/1/Data/fastq/PNA061_Sample1_S1_R1_001.fastq.gz",
      "s3://pixelgen-technologies-ngs/NextSeq2000/250314_VH00725_275_AACNHV5HV/Analysis/1/Data/fastq/PNA061_Sample2_S2_R1_001.fastq.gz"
    ),
    fastq_2 = c(
      "s3://pixelgen-technologies-ngs/NextSeq2000/250314_VH00725_275_AACNHV5HV/Analysis/1/Data/fastq/PNA061_Sample1_S1_R2_001.fastq.gz",
      "s3://pixelgen-technologies-ngs/NextSeq2000/250314_VH00725_275_AACNHV5HV/Analysis/1/Data/fastq/PNA061_Sample2_S2_R2_001.fastq.gz"
    ),
    sample_alias = c("S1", "S2"),
    condition = c("PBMC", "PHA"),
    reference_condition = c(TRUE, FALSE),
    alternative_condition = c(FALSE, TRUE)
  ), row.names = c(NA, -2L), class = c(
    "tbl_df", "tbl", "data.frame"
  ))
  es_data <- test_es_data(
    params = params,
    samplesheet = sample_sheet
  )

  expect_equal(
    default_params,
    list(
      workflow = "amplicon_demux",
      control_markers = c("mIgG1", "mIgG2a", "mIgG2b"),
      norm_method = "CLR",
      clustering_resolution = 1,
      proximity_count_cutoff = 25,
      annotation_method = "nmf",
      mc_cores = 1,
      debug_mode = FALSE,
      test_mode = FALSE
    )
  )

  expect_no_error(
    print_params(es_data)
  )

  expect_no_error(print_metadata_table(es_data))

  expect_no_error(print_session_info())

  hashed_sheet <- read_samplesheet(test_samplesheet(type = "hashing"))
  hashed_es_data <- test_es_data(
    samplesheet = hashed_sheet,
    pxl_data_processed = NULL,
    diagnostics = list(list(
      type = "pxl_load",
      target = "S1",
      message = "No PXL file was found."
    ))
  )
  hashed_table <- print_metadata_table(hashed_es_data)
  expect_equal(
    names(hashed_table$x$data),
    c("Issues", "Pool", "Sample Alias", "Sample name", "Condition")
  )
  expect_equal(
    hashed_table$x$data$Issues,
    c("\u26A0", "", "", "")
  )
})

# Create a minimal PXL-like DuckDB file with a `metadata` table and, optionally,
# a `panel_sources` table (new PXL format)
make_pxl_with_metadata <- function(metadata, panel_source_names = NULL) {
  path <- tempfile(fileext = ".pxl")
  con <- DBI::dbConnect(duckdb::duckdb(), path, read_only = FALSE)
  DBI::dbWriteTable(
    con, "metadata",
    data.frame(
      key = "sample",
      value = paste0(
        "{",
        paste(sprintf('"%s": "%s"', names(metadata), unlist(metadata)), collapse = ", "),
        "}"
      )
    )
  )
  if (!is.null(panel_source_names)) {
    n <- length(panel_source_names)
    # Rows are written in reverse so ORDER BY source_id is what fixes the order
    DBI::dbWriteTable(
      con, "panel_sources",
      data.frame(
        source_id = rev(seq_len(n) - 1L),
        name = rev(panel_source_names),
        version = "1.0.0"
      )
    )
  }
  DBI::dbDisconnect(con, shutdown = TRUE)
  path
}

test_that("print_pixelator_version works on the test data", {
  es_data <- test_es_data(
    samplesheet = read_samplesheet(test_samplesheet()),
    pxl_data_processed = get_test_data()
  )

  expect_no_error(version_table <- print_pixelator_version(es_data))
  expect_equal(
    names(version_table$x$data),
    c("Sample alias", "Sample name", "Pixelator version", "Panel name")
  )
})

test_that(".read_pxl_run_info reads the panel name from old and new PXL formats", {
  old_file <- make_pxl_with_metadata(
    list(sample_name = "S1", version = "0.30.0", panel_name = "old-panel")
  )
  new_file <- make_pxl_with_metadata(
    list(sample_name = "S2", version = "0.31.0"),
    panel_source_names = "new-panel"
  )
  multi_file <- make_pxl_with_metadata(
    list(sample_name = "S3", version = "0.31.0"),
    panel_source_names = c("panel-A", "panel-B")
  )
  missing_file <- make_pxl_with_metadata(
    list(sample_name = "S4", version = "0.31.0")
  )
  # The old format takes precedence when both are present
  both_file <- make_pxl_with_metadata(
    list(sample_name = "S5", version = "0.30.0", panel_name = "old-panel"),
    panel_source_names = "new-panel"
  )

  expect_equal(.read_pxl_run_info(old_file)$panel_name, "old-panel")
  expect_equal(.read_pxl_run_info(old_file)$version, "0.30.0")
  expect_equal(.read_pxl_run_info(new_file)$panel_name, "new-panel")
  expect_equal(.read_pxl_run_info(new_file)$sample_name, "S2")
  expect_equal(.read_pxl_run_info(multi_file)$panel_name, "panel-A+panel-B")
  expect_identical(.read_pxl_run_info(missing_file)$panel_name, NA_character_)
  expect_equal(.read_pxl_run_info(both_file)$panel_name, "old-panel")

  # Old and new files can be combined in one table
  mixed <- dplyr::bind_rows(
    list(S1 = .read_pxl_run_info(old_file), S2 = .read_pxl_run_info(new_file)),
    .id = "sample_alias"
  )
  expect_equal(mixed$panel_name, c("old-panel", "new-panel"))

  # Connections are released, so the files can be reopened for writing
  con <- DBI::dbConnect(duckdb::duckdb(), new_file, read_only = FALSE)
  expect_no_error(DBI::dbGetQuery(con, "SELECT 1"))
  DBI::dbDisconnect(con, shutdown = TRUE)
})
