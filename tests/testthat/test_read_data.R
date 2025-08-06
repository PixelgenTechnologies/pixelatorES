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

  sample_aliases <-
    c(
      "A_sample_S1" = "S1",
      "A_sample_S2" = "S2"
    )

  c(
    "run_folder/pixelator/A_sample_S1.amplicon.pxl",
    "run_folder/pixelator/A_sample_S1.demux.pxl",
    "run_folder/pixelator/A_sample_S1.graph.pxl",
    "run_folder/pixelator/A_sample_S1.denoise.pxl",
    "run_folder/pixelator/A_sample_S1.analysis.pxl",
    "run_folder/pixelator/A_sample_S1.post_analysis.pxl",
    "run_folder/pixelator/A_sample_S1.layout.pxl"
  ) %>%
    sapply(find_stage) %>%
    expect_equal(c(
      "run_folder/pixelator/A_sample_S1.amplicon.pxl" = "amplicon",
      "run_folder/pixelator/A_sample_S1.demux.pxl" = "demux",
      "run_folder/pixelator/A_sample_S1.graph.pxl" = "graph",
      "run_folder/pixelator/A_sample_S1.denoise.pxl" = "denoise",
      "run_folder/pixelator/A_sample_S1.analysis.pxl" = "analysis",
      "run_folder/pixelator/A_sample_S1.post_analysis.pxl" = "post_analysis",
      "run_folder/pixelator/A_sample_S1.layout.pxl" = "layout"
    ))

  expect_error(find_stage("run_folder/pixelator/A_sample_S1.notastage.pxl"))

  expect_no_error(res <- get_file_paths(
    file_paths = file_paths,
    sample_aliases = sample_aliases
  ))

  expect_equal(
    res,
    list(
      data_files = structure(
        list(
          sample_alias = c(
            A_sample_S1 = "S1",
            A_sample_S2 = "S2"
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
            A_sample_S1 = "S1", A_sample_S1 = "S1",
            A_sample_S1 = "S1", A_sample_S1 = "S1", A_sample_S1 = "S1",
            A_sample_S1 = "S1", A_sample_S1 = "S1", A_sample_S1 = "S1",
            A_sample_S2 = "S2", A_sample_S2 = "S2", A_sample_S2 = "S2",
            A_sample_S2 = "S2", A_sample_S2 = "S2", A_sample_S2 = "S2",
            A_sample_S2 = "S2", A_sample_S2 = "S2"
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
      )
    )
  )
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

  # Data downsampling
  expect_no_error(seur_down <- downsample_data(seur_comb,
    control_markers = c("mIgG1", "mIgG2a", "mIgG2b"),
    n_cells = 3, n_markers = 5
  ))
  expect_s4_class(seur_down, "Seurat")
  expect_equal(dim(seur_down), c(5, 6))

  # Sample sheet reading
  expect_no_error(sample_sheet <- read_samplesheet(test_samplesheet()))
  expect_equal(
    sample_sheet,
    structure(
      list(
        sample = c("S01_PBMC_unstimulated_S1", "S02_PHA_S2"),
        design = c(NA, NA),
        panel = c(NA, NA),
        fastq_1 = c(NA, NA),
        fastq_2 = c(NA, NA),
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
      sample_aliases = sample_sheet %>%
        select(sample, sample_alias) %>%
        deframe()
    ))
  expect_equal(
    data_paths$data_files %>%
      mutate(filename = str_remove(filename, ".*extdata/")),
    structure(
      list(
        sample_alias = c(
          S01_PBMC_unstimulated_S1 = "S1",
          S02_PHA_S2 = "S2"
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
          S01_PBMC_unstimulated_S1 = "S1",
          S02_PHA_S2 = "S2", S01_PBMC_unstimulated_S1 = "S1", S02_PHA_S2 = "S2",
          S01_PBMC_unstimulated_S1 = "S1", S02_PHA_S2 = "S2", S01_PBMC_unstimulated_S1 = "S1",
          S02_PHA_S2 = "S2", S01_PBMC_unstimulated_S1 = "S1", S02_PHA_S2 = "S2",
          S01_PBMC_unstimulated_S1 = "S1", S02_PHA_S2 = "S2", S01_PBMC_unstimulated_S1 = "S1",
          S02_PHA_S2 = "S2", S01_PBMC_unstimulated_S1 = "S1", S02_PHA_S2 = "S2"
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
  expect_no_error(qc_metrics <- read_qc_files(data_paths$qc_files, sample_sheet))
  expect_equal(
    names(qc_metrics),
    c("S1", "S2")
  )

  expect_equal(
    qc_metrics$S1$amplicon,
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
        vars = c("total_failed_reads", "failed_partial_upi1_umi1_reads")
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
    vars = c("q30_statistics")
  ))

  expect_no_error(
    extracted_qc_metrics <-
      extract_sample_qc_metrics(qc_metrics, "amplicon",
        vars = c(
          "a" = "total_failed_reads",
          "b" = "failed_partial_upi1_umi1_reads"
        )
      )
  )
  expect_equal(
    names(extracted_qc_metrics),
    c("sample_alias", "a", "b")
  )
})

