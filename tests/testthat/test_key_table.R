library(Seurat)

test_that("get_top_markers works as expected", {

  sample_sheet <- read_samplesheet(test_samplesheet())

  seur <-
    get_test_data()

  data_paths <-
    get_file_paths(
      data_folder = test_data_folder(),
      metadata = sample_sheet,
      sample_aliases = sample_sheet %>%
        select(sample, sample_alias) %>%
        deframe())

  sample_qc_metrics <-
    read_qc_files(data_paths$qc_files, sample_sheet)

  # Top markers
  expect_no_error(res <- get_top_markers(seur, group = "sample_alias"))
  expect_equal(
    res,
    structure(
      list(
        sample_alias = "S1",
        top3_fraction = 0.352019337952774,
        top5_fraction = 0.45234721751522,
        top_markers = "CD44, B2M, CD59, CD45, CD43"),
      row.names = c(NA, -1L),
      class = c("tbl_df", "tbl", "data.frame")
    )
  )

  # Read stats
  expect_no_error(res <- get_read_stats(seur))
  expect_equal(
    res,
    structure(
      list(sample_alias = "S1",
           n_cells = 30L,
           n_cells_over10k = 30L,
           median_reads_per_cell = 289856,
           median_abs_per_cell = 48636,
           median_isotype_count_pct = 0.0395665936206477,
           median_intracellular_count_pct = 0),
      class = c("tbl_df", "tbl", "data.frame"),
      row.names = c(NA, -1L))
  )

  # Crossing edges
  expect_no_error(res <- get_crossing_edges(sample_qc_metrics))
  expect_equal(
    res,
    structure(
      list(
        sample_alias = c("S1", "S1", "S2", "S2"),
        removed_total = c(7311081L, 7311081L, 3866617L, 3866617L),
        total_edges_in = c(144111916L, 144111916L, 90118672L, 90118672L),
        type = c("Initial stage", "Refinement stage", "Initial stage", "Refinement stage"),
        edges = c(7292268L, 18813L, 3849577L, 17040L),
        percent = c(5.06014228552759, 0.0130544374970353,
                    4.27167524173015, 0.0189084011357824)),
      row.names = c(NA, -4L),
      class = c("tbl_df", "tbl", "data.frame")
    )
  )

  # Denoising
  expect_no_error(res <- get_denoising_data(sample_qc_metrics))
  expect_equal(
    res,
    structure(
      list(
        sample_alias = c("S1", "S2"),
        ratio = c(2.67032153250738, 4.85786823694461),
        number_of_umis_removed = c(1621394L, 2724990L)),
      class = c("tbl_df", "tbl", "data.frame"),
      row.names = c(NA, -2L)
    )
  )

  # Coreness
  expect_no_error(res <- get_coreness_data(seur))
  expect_equal(
    head(res$data),
    structure(
      list(
        sample_component = c("0a45497c6bfbfb22_1", "0a45497c6bfbfb22_1",
                             "0a45497c6bfbfb22_1", "0a45497c6bfbfb22_1",
                             "0a45497c6bfbfb22_1", "2708240b908e2eba_1"),
        sample_alias = c("S1", "S1", "S1", "S1",
                         "S1", "S1"),
        molecules = c(43543L, 43543L, 43543L, 43543L, 43543L,
                      37665L),
        k_core = c("k_core_1", "k_core_2", "k_core_3", "k_core_4",
                   "k_core_5", "k_core_1"),
        nodes = c(10260, 6046, 19434, 7803,
                  0, 9596),
        percent_nodes = c(23.5629148198333, 13.8851250488023,
                          44.6317433341754, 17.920216797189,
                          0, 25.4772335059073),
        coreness = c(1,
                     2, 3, 4, 5, 1)),
      row.names = c(NA, -6L),
      class = c("tbl_df", "tbl", "data.frame")
    )
  )

  expect_equal(
    head(res$component_summary),
    structure(
      list(
        sample_component = c("0a45497c6bfbfb22_1", "0a45497c6bfbfb22_2",
                             "0a45497c6bfbfb22_3", "0a45497c6bfbfb22_4",
                             "0a45497c6bfbfb22_5", "0a45497c6bfbfb22_6"),
        sample_alias = c("S1", "S1", "S1", "S1",
                         "S1", "S1"),
        mean_coreness = c(2.5690926210872, 2.5690926210872,
                          2.5690926210872, 2.5690926210872,
                          2.5690926210872, 2.5690926210872
        ),
        percent_dangling_nodes = c(23.5629148198333, 23.5629148198333,
                                   23.5629148198333, 23.5629148198333,
                                   23.5629148198333, 23.5629148198333
        ),
        percent_well_connected_nodes = c(62.5519601313644, 62.5519601313644,
                                         62.5519601313644, 62.5519601313644,
                                         62.5519601313644, 62.5519601313644
        )),
      class = c("grouped_df", "tbl_df", "tbl", "data.frame"),
      row.names = c(NA, -6L),
      groups = structure(
        list(
          sample_component = c("0a45497c6bfbfb22_1", "0a45497c6bfbfb22_2",
                               "0a45497c6bfbfb22_3", "0a45497c6bfbfb22_4",
                               "0a45497c6bfbfb22_5", "0a45497c6bfbfb22_6"),
          .rows = structure(list(1L, 2L, 3L, 4L, 5L, 6L),
                            ptype = integer(0),
                            class = c("vctrs_list_of", "vctrs_vctr", "list"))),
        class = c("tbl_df", "tbl", "data.frame"),
        row.names = c(NA, -6L),
        .drop = TRUE))
  )

  expect_equal(
    res$sample_summary,
    structure(
      list(
        sample_alias = "S1",
        median_mean_coreness = 2.47440345751157,
        median_percent_dangling_nodes = 24.993831729583,
        median_percent_well_connected_nodes = 57.1067689310933),
      class = c("tbl_df", "tbl", "data.frame"),
      row.names = c(NA, -1L)
    )
  )

})
