test_that("File location works as expected", {
  file_paths <-
    c(
      "run_folder/pixelator/amplicon/A_sample_S1.report.json",
      "run_folder/pixelator/analysis/A_sample_S1.report.json",
      "run_folder/pixelator/collapse/A_sample_S1.collapse.m1.part_000.report.json",
      "run_folder/pixelator/collapse/A_sample_S1.collapse.m1.part_001.report.json",
      "run_folder/pixelator/collapse/A_sample_S1.collapse.m2.part_000.report.json",
      "run_folder/pixelator/collapse/A_sample_S1.collapse.m2.part_001.report.json",
      "run_folder/pixelator/collapse/A_sample_S1.report.json",
      "run_folder/pixelator/demux/A_sample_S1.report.json",
      "run_folder/pixelator/denoise/A_sample_S1.report.json",
      "run_folder/pixelator/graph/A_sample_S1.report.json",
      "run_folder/pixelator/layout/A_sample_S1.report.json",
      "run_folder/pixelator/A_sample_S1.layout.pxl",
      "run_folder/pixelator/post_analysis/A_sample_S1.report.json",
      "run_folder/pixelator/amplicon/A_sample_S2.report.json",
      "run_folder/pixelator/analysis/A_sample_S2.report.json",
      "run_folder/pixelator/collapse/A_sample_S2.collapse.m1.part_000.report.json",
      "run_folder/pixelator/collapse/A_sample_S2.collapse.m1.part_001.report.json",
      "run_folder/pixelator/collapse/A_sample_S2.collapse.m2.part_000.report.json",
      "run_folder/pixelator/collapse/A_sample_S2.collapse.m2.part_001.report.json",
      "run_folder/pixelator/collapse/A_sample_S2.report.json",
      "run_folder/pixelator/demux/A_sample_S2.report.json",
      "run_folder/pixelator/denoise/A_sample_S2.report.json",
      "run_folder/pixelator/graph/A_sample_S2.report.json",
      "run_folder/pixelator/layout/A_sample_S2.report.json",
      "run_folder/pixelator/A_sample_S2.layout.pxl",
      "run_folder/pixelator/post_analysis/A_sample_S2.report.json"
    )

  sample_sheet <-
    structure(list(sample = c("A_sample_S1", "A_sample_S2"), sample_alias = c("S1", "S2"), condition = c("PBMC", "PHA")), row.names = c(
      NA,
      -2L
    ), class = c("tbl_df", "tbl", "data.frame"))

  c(
    "run_folder/pixelator/A_sample_S1.amplicon.pxl",
    "run_folder/pixelator/A_sample_S1.demux.pxl",
    "run_folder/pixelator/A_sample_S1.graph.pxl",
    "run_folder/pixelator/A_sample_S1.denoise.pxl",
    "run_folder/pixelator/A_sample_S1.analysis.pxl",
    "run_folder/pixelator/A_sample_S1.post_analysis.pxl",
    "run_folder/pixelator/A_sample_S1.layout.pxl"
  ) %>%
    sapply(find_stage, stages = amplicon_stages()$all) %>%
    expect_equal(c(
      "run_folder/pixelator/A_sample_S1.amplicon.pxl" = "amplicon",
      "run_folder/pixelator/A_sample_S1.demux.pxl" = "demux",
      "run_folder/pixelator/A_sample_S1.graph.pxl" = "graph",
      "run_folder/pixelator/A_sample_S1.denoise.pxl" = "denoise",
      "run_folder/pixelator/A_sample_S1.analysis.pxl" = "analysis",
      "run_folder/pixelator/A_sample_S1.post_analysis.pxl" = "post_analysis",
      "run_folder/pixelator/A_sample_S1.layout.pxl" = "layout"
    ))

  expect_error(find_stage(
    "run_folder/pixelator/A_sample_S1.notastage.pxl",
    stages = amplicon_stages()$all
  ))
  expect_error(find_stage("run_folder/pixelator/A_sample_S1.graph.pxl"))

  expect_no_error(res <- get_file_paths(
    file_paths = file_paths,
    sample_sheet = sample_sheet,
    stages = amplicon_stages()
  ))

  expect_equal(
    res,
    list(
      data_files = structure(
        list(
          sample_alias = c(
            "S1",
            "S2"
          ),
          filename = c(
            "run_folder/pixelator/A_sample_S1.layout.pxl",
            "run_folder/pixelator/A_sample_S2.layout.pxl"
          )
        ),
        row.names = c(
          NA,
          -2L
        ),
        class = c("tbl_df", "tbl", "data.frame")
      ),
      qc_files = structure(
        list(
          sample_alias = c(
            "S1", "S1",
            "S1", "S1", "S1",
            "S1", "S1", "S1",
            "S2", "S2", "S2",
            "S2", "S2", "S2",
            "S2", "S2"
          ),
          filename = c(
            "run_folder/pixelator/amplicon/A_sample_S1.report.json",
            "run_folder/pixelator/analysis/A_sample_S1.report.json",
            "run_folder/pixelator/collapse/A_sample_S1.report.json",
            "run_folder/pixelator/demux/A_sample_S1.report.json",
            "run_folder/pixelator/denoise/A_sample_S1.report.json",
            "run_folder/pixelator/graph/A_sample_S1.report.json",
            "run_folder/pixelator/layout/A_sample_S1.report.json",
            "run_folder/pixelator/post_analysis/A_sample_S1.report.json",
            "run_folder/pixelator/amplicon/A_sample_S2.report.json",
            "run_folder/pixelator/analysis/A_sample_S2.report.json",
            "run_folder/pixelator/collapse/A_sample_S2.report.json",
            "run_folder/pixelator/demux/A_sample_S2.report.json",
            "run_folder/pixelator/denoise/A_sample_S2.report.json",
            "run_folder/pixelator/graph/A_sample_S2.report.json",
            "run_folder/pixelator/layout/A_sample_S2.report.json",
            "run_folder/pixelator/post_analysis/A_sample_S2.report.json"
          ),
          stage = c(
            `run_folder/pixelator/amplicon/A_sample_S1.report.json` =
              "amplicon",
            `run_folder/pixelator/analysis/A_sample_S1.report.json` =
              "analysis",
            `run_folder/pixelator/collapse/A_sample_S1.report.json` =
              "collapse",
            `run_folder/pixelator/demux/A_sample_S1.report.json` =
              "demux",
            `run_folder/pixelator/denoise/A_sample_S1.report.json` =
              "denoise",
            `run_folder/pixelator/graph/A_sample_S1.report.json` =
              "graph",
            `run_folder/pixelator/layout/A_sample_S1.report.json` =
              "layout",
            `run_folder/pixelator/post_analysis/A_sample_S1.report.json` =
              "post_analysis",
            `run_folder/pixelator/amplicon/A_sample_S2.report.json` =
              "amplicon",
            `run_folder/pixelator/analysis/A_sample_S2.report.json` =
              "analysis",
            `run_folder/pixelator/collapse/A_sample_S2.report.json` =
              "collapse",
            `run_folder/pixelator/demux/A_sample_S2.report.json` =
              "demux",
            `run_folder/pixelator/denoise/A_sample_S2.report.json` =
              "denoise",
            `run_folder/pixelator/graph/A_sample_S2.report.json` =
              "graph",
            `run_folder/pixelator/layout/A_sample_S2.report.json` =
              "layout",
            `run_folder/pixelator/post_analysis/A_sample_S2.report.json` =
              "post_analysis"
          )
        ),
        row.names = c(NA, -16L), class = c("tbl_df", "tbl", "data.frame")
      ),
      pool_qc_files = NULL
    )
  )
})

test_that("Custom stage vocabularies drive file discovery as expected", {
  custom_stages <- list(
    all = c("amplicon", "demux", "graph", "analysis"),
    pool = "amplicon",
    pxl_preference = c("analysis", "graph")
  )

  expect_equal(
    find_stage("run/A_sample_S1.graph.pxl", stages = custom_stages$all),
    "graph"
  )
  expect_error(
    find_stage("run/A_sample_S1.layout.pxl", stages = custom_stages$all)
  )
  expect_equal(
    find_stage("run/A_sample_S1.layout.pxl", stages = amplicon_stages()$all),
    "layout"
  )

  file_paths <- c(
    "run/amplicon/A_sample_S1.report.json",
    "run/graph/A_sample_S1.report.json",
    "run/analysis/A_sample_S1.report.json",
    "run/A_sample_S1.graph.pxl",
    "run/A_sample_S1.analysis.pxl"
  )
  sample_sheet <- tibble(
    sample = "A_sample_S1",
    sample_alias = "S1",
    condition = "PBMC"
  )

  # pxl_preference decides the winner: the custom vocabulary ranks graph last.
  expect_equal(
    get_file_paths(
      file_paths = file_paths,
      sample_sheet = sample_sheet,
      stages = custom_stages
    )$data_files$filename,
    "run/A_sample_S1.graph.pxl"
  )
  expect_equal(
    get_file_paths(
      file_paths = file_paths,
      sample_sheet = sample_sheet,
      stages = amplicon_stages()
    )$data_files$filename,
    "run/A_sample_S1.analysis.pxl"
  )

  expect_error(
    get_file_paths(
      file_paths = file_paths,
      sample_sheet = sample_sheet,
      stages = list(all = "graph")
    )
  )

  # The collapse shard filter only applies to vocabularies that have the stage.
  shard_paths <- c(
    "run/collapse/A_sample_S1.collapse.m1.part_000.report.json",
    "run/collapse/A_sample_S1.report.json",
    "run/graph/A_sample_S1.report.json",
    "run/A_sample_S1.graph.pxl"
  )

  expect_equal(
    nrow(get_file_paths(
      file_paths = shard_paths,
      sample_sheet = sample_sheet,
      stages = amplicon_stages()
    )$qc_files),
    2
  )
  expect_equal(
    nrow(get_file_paths(
      file_paths = shard_paths,
      sample_sheet = sample_sheet,
      allow_unknown = TRUE,
      stages = custom_stages
    )$qc_files),
    3
  )
})

test_that("Custom stage vocabularies drive QC extraction as expected", {
  qc_metrics <- list(
    qc_files = list(
      S1 = list(build = list(n = 1L), tally = list(n = 2L)),
      S2 = list(build = list(n = 3L), tally = list(n = 4L))
    ),
    pool_qc_files = list(
      P1 = list(build = list(n = 98L), tally = list(n = 99L))
    )
  )
  custom_stages <- list(
    all = c("build", "tally"),
    pool = "tally",
    pxl_preference = "build"
  )

  expect_equal(
    extract_sample_qc_metrics(
      qc_metrics,
      vars = "n",
      stage = "build",
      stages = custom_stages
    ),
    tibble(sample_alias = c("S1", "S2"), n = c(1L, 3L))
  )

  # `pool` routes the stage to pool-level QC instead.
  expect_equal(
    extract_sample_qc_metrics(
      qc_metrics,
      vars = "n",
      stage = "tally",
      stages = custom_stages
    ),
    tibble(sample_alias = "P1", n = 99L)
  )

  expect_error(
    extract_sample_qc_metrics(
      qc_metrics,
      vars = "n",
      stage = "graph",
      stages = custom_stages
    )
  )
  expect_error(
    extract_sample_qc_metrics(
      qc_metrics,
      vars = "n",
      stage = "build",
      stages = amplicon_stages()
    )
  )
  expect_error(extract_sample_qc_metrics(qc_metrics, vars = "n", stage = "build"))
})


test_that("File reading works as expected", {
  data_files <-
    structure(
      list(
        sample_alias = c(
          S01_PBMC_unstimulated_S1 = "S1",
          S02_PHA_S2 = "S2"
        ),
        filename = c(
          minimal_pna_pxl_file(),
          minimal_pna_pxl_file()
        )
      ),
      row.names = c(NA, -2L), class = c("tbl_df", "tbl", "data.frame")
    )

  data_folder <-
    data_files$filename[1] %>%
    dirname()

  sample_sheet <-
    structure(list(
      sample = c("S01_PBMC_unstimulated_S1", "S02_PHA_S2"), design = c("pna-2", "pna-2"), panel = c(
        "pna-rnd-158plex-final",
        "pna-rnd-158plex-final"
      ), fastq_1 = c(
        "s3://fastq/PNA061_Sample1_S1_R1_001.fastq.gz",
        "s3://fastq/PNA061_Sample2_S2_R1_001.fastq.gz"
      ), fastq_2 = c(
        "s3://fastq/PNA061_Sample1_S1_R2_001.fastq.gz",
        "s3://fastq/PNA061_Sample2_S2_R2_001.fastq.gz"
      ), sample_alias = c("S1", "S2"), condition = c("PBMC", "PHA")
    ), row.names = c(NA, -2L), class = c(
      "tbl_df", "tbl",
      "data.frame"
    ))

  expect_no_error(seur_list <- load_pxl_data_list(data_folder, data_files, sample_sheet))
  expect_s4_class(seur_list[[1]], "Seurat")

  expect_error(load_pxl_data_list(data_files, sample_sheet[1, ]))
  expect_error(load_pxl_data_list(data_files[1, ], sample_sheet))

  # Data merging
  expect_no_error(seur_comb <- merge_data(seur_list, sample_sheet))
  expect_s4_class(seur_comb, "Seurat")
  expect_equal(dim(seur_comb), dim(seur_list[[1]]) * c(1, 2))

  expect_error(merge_data(unname(seur_list), sample_sheet))
  expect_error(merge_data(unname(seur_list[1]), sample_sheet))

  # Data downsampling
  expect_no_error(seur_down <- downsample_data(seur_comb,
    control_markers = c("mIgG1", "mIgG2a", "mIgG2b"),
    n_cells = 3, n_markers = 5
  ))
  expect_s4_class(seur_down, "Seurat")
  expect_equal(dim(seur_down), c(5, 6))

  # Edge case 1: fewer cells than requested in at least one sample
  expect_warning(
    seur_down_low_cells <- downsample_data(
      seur_comb,
      control_markers = c("mIgG1", "mIgG2a", "mIgG2b"),
      n_cells = 1000,
      n_markers = 5
    ),
    "fewer cells"
  )
  expect_equal(ncol(seur_down_low_cells), ncol(seur_comb))

  # Edge case 2: fewer non-control markers available than requested
  all_markers <- rownames(seur_comb)
  expect_gt(length(all_markers), 1)
  control_set <- all_markers[-length(all_markers)]
  expect_warning(
    seur_down_low_markers <- downsample_data(
      seur_comb,
      control_markers = control_set,
      n_cells = 3,
      n_markers = nrow(seur_comb) + 5
    ),
    "non-control markers"
  )
  expect_equal(nrow(seur_down_low_markers), nrow(seur_comb))

  # Sample sheet reading
  expect_no_error(sample_sheet <- read_samplesheet(test_samplesheet()))
  expect_equal(
    sample_sheet,
    structure(
      list(
        sample = c("S01_PBMC_unstimulated_S1", "S02_PHA_S2"),
        sample_alias = c("S1", "S2"),
        condition = c("PBMC", "PHA")
      ),
      row.names = c(NA, -2L),
      class = c(
        "tbl_df", "tbl",
        "data.frame"
      )
    )
  )


  expect_no_error(data_paths <-
    get_file_paths(
      data_folder = test_data_folder(),
      sample_sheet = sample_sheet,
      stages = amplicon_stages()
    ))
  expect_equal(
    data_paths$data_files %>%
      mutate(filename = str_remove(filename, ".*extdata/")),
    structure(
      list(
        sample_alias = c(
          "S1",
          "S2"
        ),
        filename = c(
          "qc_jsons/S01_PBMC_unstimulated_S1.layout.pxl",
          "qc_jsons/S02_PHA_S2.layout.pxl"
        )
      ),
      row.names = c(NA, -2L),
      class = c(
        "tbl_df",
        "tbl", "data.frame"
      )
    )
  )
  expect_equal(
    data_paths$qc_files %>%
      mutate(
        filename = str_remove(filename, ".*extdata/"),
        stage = unname(stage)
      ),
    structure(
      list(
        sample_alias = c(
          "S1",
          "S2", "S1", "S2",
          "S1", "S2", "S1",
          "S2", "S1", "S2",
          "S1", "S2", "S1",
          "S2", "S1", "S2"
        ), filename = c(
          "qc_jsons/amplicon/S01_PBMC_unstimulated_S1.report.json",
          "qc_jsons/amplicon/S02_PHA_S2.report.json", "qc_jsons/analysis/S01_PBMC_unstimulated_S1.report.json",
          "qc_jsons/analysis/S02_PHA_S2.report.json", "qc_jsons/collapse/S01_PBMC_unstimulated_S1.report.json",
          "qc_jsons/collapse/S02_PHA_S2.report.json", "qc_jsons/demux/S01_PBMC_unstimulated_S1.report.json",
          "qc_jsons/demux/S02_PHA_S2.report.json", "qc_jsons/denoise/S01_PBMC_unstimulated_S1.report.json",
          "qc_jsons/denoise/S02_PHA_S2.report.json", "qc_jsons/graph/S01_PBMC_unstimulated_S1.report.json",
          "qc_jsons/graph/S02_PHA_S2.report.json", "qc_jsons/layout/S01_PBMC_unstimulated_S1.report.json",
          "qc_jsons/layout/S02_PHA_S2.report.json", "qc_jsons/post_analysis/S01_PBMC_unstimulated_S1.report.json",
          "qc_jsons/post_analysis/S02_PHA_S2.report.json"
        ),
        stage = c(
          "amplicon",
          "amplicon", "analysis", "analysis", "collapse", "collapse", "demux",
          "demux", "denoise", "denoise", "graph", "graph", "layout", "layout",
          "post_analysis", "post_analysis"
        )
      ),
      row.names = c(NA, -16L),
      class = c(
        "tbl_df",
        "tbl", "data.frame"
      )
    )
  )

  # QC data reading
  expect_no_error(qc_metrics <- read_qc_files(data_paths, sample_sheet))
  expect_equal(
    names(qc_metrics),
    c("qc_files", "pool_qc_files")
  )

  expect_equal(
    qc_metrics$qc_files$S1$amplicon,
    list(
      sample_id = "S01_PBMC_unstimulated_S1", product_id = "single-cell-pna",
      report_type = "amplicon", input_reads = 400745924L, output_reads = 365893238L,
      passed_missing_uei_reads = 54042998L, passed_partial_uei_reads = 6355010L,
      passed_missing_lbs1_anchor = 31309030L, failed_too_many_n_reads = 1077746L,
      failed_partial_upi1_umi1_reads = 0L, failed_partial_upi2_umi2_reads = 0L,
      failed_missing_upi1_umi1_reads = 22721440L, failed_missing_upi2_umi2_reads = 11053500L,
      total_failed_reads = 34852686L,
      q30_statistics = list(
        total = 0.957683900982373,
        umi1 = 0.972854977889167, pid1 = 0.976935098483564, lbs1 = 0.943076352015842,
        uei = 0.827510713020246, lbs2 = 0.884484284209117, pid2 = 0.970777786825347,
        umi2 = 0.97077052793588
      ),
      basepair_counts = list(
        input = 48891002728,
        input_read1 = 17632820656, input_read2 = 31258182072,
        quality_trimmed = 505928217L, quality_trimmed_read1 = 392817490L,
        quality_trimmed_read2 = 113110727L, output = 51956839796
      ),
      failed_invalid_amplicon_reads = 33774940L, fraction_discarded_reads = 0.0869695333445238
    )
  )

  expect_no_error(
    extracted_qc_metrics <-
      extract_sample_qc_metrics(qc_metrics, "amplicon",
        vars = c("total_failed_reads", "failed_partial_upi1_umi1_reads"),
        stages = amplicon_stages()
      )
  )

  expect_equal(
    extracted_qc_metrics,
    structure(
      list(
        sample_alias = c("S1", "S2"),
        total_failed_reads = c(
          34852686L,
          18770743L
        ),
        failed_partial_upi1_umi1_reads = c(0L, 0L)
      ),
      class = c(
        "tbl_df",
        "tbl", "data.frame"
      ),
      row.names = c(NA, -2L)
    )
  )
  expect_error(extract_sample_qc_metrics(qc_metrics, "amplicon",
    vars = c("q30_statistics"),
    stages = amplicon_stages()
  ))

  expect_no_error(
    extracted_qc_metrics <-
      extract_sample_qc_metrics(qc_metrics, "amplicon",
        vars = c(
          "a" = "total_failed_reads",
          "b" = "failed_partial_upi1_umi1_reads"
        ),
        stages = amplicon_stages()
      )
  )
  expect_equal(
    names(extracted_qc_metrics),
    c("sample_alias", "a", "b")
  )


  # Hashing
  expect_no_error(sample_sheet_hashing <- read_samplesheet(test_samplesheet(type = "hashing")))

  expect_equal(
    sample_sheet_hashing,
    structure(
      list(
        sample = c("S1", "S2", "S11", "S12"),
        sample_alias = c("S1", "S2", "S11", "S12"),
        condition = c("PBMC", "Raji", "PBMC", "Raji"),
        pool = c("pool1", "pool1", "pool2", "pool2")
      ),
      row.names = c(NA, -4L),
      class = c("tbl_df", "tbl", "data.frame")
    )
  )

  # Error out when duplicates between pool and samples
  temp_file <- tempfile(fileext = ".csv")
  sample_sheet_hashing %>%
    mutate(pool = sample) %>%
    write.csv(temp_file)
  expect_error(
    suppressMessages(read_samplesheet(temp_file)),
    regexp = "The `pool` column.*cannot contain.*also present"
  )

  expect_no_error(data_paths <-
    get_file_paths(
      data_folder = test_data_folder(type = "hashing"),
      sample_sheet = sample_sheet_hashing,
      stages = amplicon_stages()
    ))

  expect_equal(
    data_paths$data_files %>%
      mutate(filename = str_remove(filename, ".*extdata/")),
    structure(list(sample_alias = c("S1", "S11", "S12", "S2"), filename = c(
      "qc_jsons_hashing/S1.layout.pxl",
      "qc_jsons_hashing/S11.layout.pxl", "qc_jsons_hashing/S12.layout.pxl",
      "qc_jsons_hashing/S2.layout.pxl"
    )), row.names = c(NA, -4L), class = c(
      "tbl_df",
      "tbl", "data.frame"
    ))
  )
  expect_equal(
    data_paths$qc_files %>%
      mutate(
        filename = str_remove(filename, ".*extdata/"),
        stage = unname(stage)
      ),
    structure(list(sample_alias = c(
      "S1", "S11", "S12", "S2", "S1",
      "S11", "S12", "S2", "S1", "S11", "S12", "S2", "S1", "S11", "S12",
      "S2", "S1", "S11", "S12", "S2"
    ), filename = c(
      "qc_jsons_hashing/analysis/S1.report.json",
      "qc_jsons_hashing/analysis/S11.report.json", "qc_jsons_hashing/analysis/S12.report.json",
      "qc_jsons_hashing/analysis/S2.report.json", "qc_jsons_hashing/denoise/S1.report.json",
      "qc_jsons_hashing/denoise/S11.report.json", "qc_jsons_hashing/denoise/S12.report.json",
      "qc_jsons_hashing/denoise/S2.report.json", "qc_jsons_hashing/layout/S1.report.json",
      "qc_jsons_hashing/layout/S11.report.json", "qc_jsons_hashing/layout/S12.report.json",
      "qc_jsons_hashing/layout/S2.report.json", "qc_jsons_hashing/post_analysis/S1.report.json",
      "qc_jsons_hashing/post_analysis/S11.report.json", "qc_jsons_hashing/post_analysis/S12.report.json",
      "qc_jsons_hashing/post_analysis/S2.report.json", "qc_jsons_hashing/sample_calling/S1.report.json",
      "qc_jsons_hashing/sample_calling/S11.report.json", "qc_jsons_hashing/sample_calling/S12.report.json",
      "qc_jsons_hashing/sample_calling/S2.report.json"
    ), stage = c(
      "analysis",
      "analysis", "analysis", "analysis", "denoise", "denoise", "denoise",
      "denoise", "layout", "layout", "layout", "layout", "post_analysis",
      "post_analysis", "post_analysis", "post_analysis", "sample_calling",
      "sample_calling", "sample_calling", "sample_calling"
    )), row.names = c(
      NA,
      -20L
    ), class = c("tbl_df", "tbl", "data.frame"))
  )

  expect_error(
    read_qc_files(
      data_paths,
      structure(
        list(
          sample = c("S01_PBMC_unstimulated_S1", "S02_PHA_S2", "S3"),
          sample_alias = c("S1", "S2", "S3"), condition = c("PBMC", "PHA", "S3")
        ),
        row.names = c(NA, -3L), class = c("tbl_df", "tbl", "data.frame")
      )
    )
  )
  # QC data reading
  expect_no_error(qc_metrics <- read_qc_files(data_paths, sample_sheet_hashing))
  expect_equal(
    names(qc_metrics),
    c("qc_files", "pool_qc_files")
  )

  expect_equal(
    qc_metrics$pool_qc_files$pool1$amplicon,
    list(
      sample_id = "pool1", product_id = "single-cell-pna", report_type = "amplicon",
      input_reads = 2240712L, output_reads = 1974950L, passed_missing_uei_reads = 102163L,
      passed_partial_uei_reads = 41185L, passed_missing_lbs1_anchor = 37884L,
      failed_too_many_n_reads = 70306L, failed_partial_upi1_umi1_reads = 0L,
      failed_partial_upi2_umi2_reads = 0L, failed_missing_upi1_umi1_reads = 124069L,
      failed_missing_upi2_umi2_reads = 69176L, failed_lbs_detected_in_umi_reads = 1868L,
      failed_low_complexity_umi_reads = 343L, total_failed_reads = 265762L,
      q30_statistics = list(
        total = 0.96389646747869, umi1 = 0.97554708437465,
        pid1 = 0.978029064026937, lbs1 = 0.950325084457053, uei = 0.915778762331536,
        lbs2 = 0.926549302757264, pid2 = 0.969809109091369, umi2 = 0.972393044308536
      ),
      basepair_counts = list(
        input = 273366864L, input_read1 = 98591328L,
        input_read2 = 174775536L, quality_trimmed = 2770353L,
        quality_trimmed_read1 = 2168472L, quality_trimmed_read2 = 601881L,
        output = 280442900L
      ), failed_invalid_amplicon_reads = 193245L,
      fraction_discarded_reads = 0.118606050219752
    )
  )

  expect_no_error(
    extracted_qc_metrics <-
      extract_sample_qc_metrics(qc_metrics,
        stage = "amplicon",
        vars = c("total_failed_reads", "failed_partial_upi1_umi1_reads"),
        stages = amplicon_stages()
      )
  )

  expect_equal(
    extracted_qc_metrics,
    structure(list(sample_alias = c("pool1", "pool2"), total_failed_reads = c(
      265762L,
      265762L
    ), failed_partial_upi1_umi1_reads = c(0L, 0L)), row.names = c(
      NA,
      -2L
    ), class = c("tbl_df", "tbl", "data.frame"))
  )
  expect_error(extract_sample_qc_metrics(qc_metrics, "amplicon",
    vars = c("q30_statistics"),
    stages = amplicon_stages()
  ))

  expect_no_error(
    extracted_qc_metrics <-
      extract_sample_qc_metrics(qc_metrics, "amplicon",
        vars = c(
          "a" = "total_failed_reads",
          "b" = "failed_partial_upi1_umi1_reads"
        ),
        stages = amplicon_stages()
      )
  )
  expect_equal(
    names(extracted_qc_metrics),
    c("sample_alias", "a", "b")
  )
})
