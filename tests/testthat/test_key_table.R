library(Seurat)

data_types <- c("default", "hashing")

for (data_type in data_types) {
  test_that(paste("get_top_markers works as expected for", data_type, "data"), {
    sample_sheet <- read_samplesheet(test_samplesheet(type = data_type))

    seur <-
      get_test_data(type = data_type)

    data_paths <-
      get_file_paths(
        data_folder = test_data_folder(type = data_type),
        sample_sheet = sample_sheet
      )

    sample_qc_metrics <-
      read_qc_files(data_paths, sample_sheet)

    # Top markers
    expect_no_error(res <- get_top_markers(seur, group = "sample_alias"))
    expect_equal(
      res,
      structure(
        list(
          sample_alias = "S1",
          top3_fraction = 0.352019337952774,
          top5_fraction = 0.45234721751522,
          top_markers = "CD44, B2M, CD59, CD45, CD43"
        ),
        row.names = c(NA, -1L),
        class = c("tbl_df", "tbl", "data.frame")
      )
    )

    # Sequencing saturation
    expect_no_error(res <- get_seq_saturation(seur, sample_qc_metrics))
    expect_equal(
      res,
      switch(data_type,
             default =
               structure(list(
                 sample_alias = "S1", Q30 = 0.957683900982373,
                 total_reads = 400745924L, deduped_valid_reads = 232530522L,
                 valid_reads = 359952911L, graph_edges = 3171564L, graph_proteins = 1469442L,
                 graph_reads = 8746572L, valid_reads_saturation = 35.3997384396761,
                 graph_edge_saturation = 63.7393483984354, graph_node_saturation = 83.1997953026626,
                 fraction_valid_reads = 89.8207291560625, fraction_graph_reads = 2.18257291620014
               ), row.names = c(
                 NA,
                 -1L
               ), class = c("tbl_df", "tbl", "data.frame")),
             hashing =
               structure(list(
                 pool = "pool1", total_reads = 2240712L, deduped_valid_reads = 1159113L,
                 valid_reads = 1953716L, graph_edges = 3171564L, graph_proteins = 1469442L,
                 graph_reads = 8746572L, valid_reads_saturation = 40.6713667697864,
                 graph_edge_saturation = 63.7393483984354, graph_node_saturation = 83.1997953026626,
                 fraction_valid_reads = 87.1917497652532, fraction_graph_reads = 390.347889420863
               ), row.names = c(
                 NA,
                 -1L
               ), class = c("tbl_df", "tbl", "data.frame"))
      )
    )

    # Read stats
    expect_no_error(res <- get_read_stats(seur))
    expect_equal(
      res,
      structure(
        list(
          sample_alias = "S1",
          n_cells = 30L,
          n_cells_over10k = 30L,
          median_reads_per_cell = 289856,
          median_abs_per_cell = 48636,
          median_isotype_count_pct = 0.0395665936206477,
          median_intracellular_count_pct = 0
        ),
        class = c("tbl_df", "tbl", "data.frame"),
        row.names = c(NA, -1L)
      )
    )

    # Crossing edges
    expect_no_error(res <- get_crossing_edges(sample_qc_metrics))
    expect_equal(
      res,
      switch(data_type,
             default =
               structure(
                 list(
                   sample_alias = c("S1", "S1", "S2", "S2"),
                   removed_total = c(7311081L, 7311081L, 3866617L, 3866617L),
                   total_edges_in = c(144111916L, 144111916L, 90118672L, 90118672L),
                   type = c("Initial stage", "Refinement stage", "Initial stage", "Refinement stage"),
                   edges = c(7292268L, 18813L, 3849577L, 17040L),
                   percent = c(
                     5.06014228552759, 0.0130544374970353,
                     4.27167524173015, 0.0189084011357824
                   )
                 ),
                 row.names = c(NA, -4L),
                 class = c("tbl_df", "tbl", "data.frame")
               ),
             hashing =
               structure(list(
                 pool = c("pool1", "pool1", "pool2", "pool2"),
                 removed_total = c(13678L, 13678L, 13564L, 13564L), total_edges_in = c(
                   813120L,
                   813120L, 813117L, 813117L
                 ), type = c(
                   "Initial stage", "Refinement stage",
                   "Initial stage", "Refinement stage"
                 ), edges = c(
                   9910L, 3768L,
                   10008L, 3556L
                 ), percent = c(
                   1.21876229830775, 0.463400236127509,
                   1.23081918100347, 0.437329437215063
                 )
               ), row.names = c(
                 NA,
                 -4L
               ), class = c("tbl_df", "tbl", "data.frame"))
      )
    )

    # Degree distribution

    expect_no_error(res <- get_degree_distribution(sample_qc_metrics))
    expect_equal(
      head(res),
      switch (data_type,
              default =
                structure(list(sample_alias = c("S1", "S1", "S1", "S1", "S1",
                                                "S1"), umi_type = c("umi1", "umi1", "umi1", "umi1", "umi1", "umi1"
                                                ), degree = 1:6, n = c(27439568L, 7545696L, 5571202L, 4951939L,
                                                                       4609669L, 4290856L), percent_nodes = c(19.0277655876346, 5.23250711102857,
                                                                                                              3.86330884281272, 3.43388549325068, 3.19654089191878, 2.97546237383531
                                                                       )), row.names = c(NA, -6L), class = c("tbl_df", "tbl", "data.frame"
                                                                       )),
              hashing =
                structure(list(pool = c("pool1", "pool1", "pool1", "pool1", "pool1",
                                        "pool1"), umi_type = c("umi1", "umi1", "umi1", "umi1", "umi1",
                                                               "umi1"), degree = 1:6, n = c(75201L, 55220L, 45622L, 33576L,
                                                                                            22807L, 14510L), percent_nodes = c(14.6206175184554, 10.7359011099467,
                                                                                                                               8.86985295976078, 6.52786337681224, 4.43414879780071, 2.82104174402983
                                                                                            )), row.names = c(NA, -6L), class = c("tbl_df", "tbl", "data.frame"
                                                                                            ))
    )
    )

    # Denoising
    expect_no_error(res <- get_denoising_data(sample_qc_metrics))
    expect_equal(
      res,
      switch(data_type,
             default =
               structure(
                 list(
                   sample_alias = c("S1", "S2"),
                   ratio = c(2.67032153250738, 4.85786823694461),
                   number_of_umis_removed = c(1621394L, 2724990L)
                 ),
                 class = c("tbl_df", "tbl", "data.frame"),
                 row.names = c(NA, -2L)
               ),
             hashing =
               structure(list(sample_alias = c("S1", "S11", "S12", "S2"), ratio = c(
                 4.78266613819527,
                 4.78980422103143, 4.73014521545811, 4.6124845823863
               ), number_of_umis_removed = c(
                 3620L,
                 3638L, 6000L, 5946L
               )), row.names = c(NA, -4L), class = c(
                 "tbl_df",
                 "tbl", "data.frame"
               ))
      )
    )

    # Coreness
    expect_no_error(res <- get_coreness_data(seur))
    expect_equal(
      head(res$data),
      structure(
        list(
          sample_component = c(
            "0a45497c6bfbfb22_1", "0a45497c6bfbfb22_1",
            "0a45497c6bfbfb22_1", "0a45497c6bfbfb22_1",
            "0a45497c6bfbfb22_1", "2708240b908e2eba_1"
          ),
          sample_alias = c(
            "S1", "S1", "S1", "S1",
            "S1", "S1"
          ),
          molecules = c(
            43543L, 43543L, 43543L, 43543L, 43543L,
            37665L
          ),
          k_core = c(
            "k_core_1", "k_core_2", "k_core_3", "k_core_4",
            "k_core_5", "k_core_1"
          ),
          nodes = c(
            10260, 6046, 19434, 7803,
            0, 9596
          ),
          percent_nodes = c(
            23.5629148198333, 13.8851250488023,
            44.6317433341754, 17.920216797189,
            0, 25.4772335059073
          ),
          coreness = c(
            1,
            2, 3, 4, 5, 1
          )
        ),
        row.names = c(NA, -6L),
        class = c("tbl_df", "tbl", "data.frame")
      )
    )

    expect_equal(
      head(res$component_summary),
      structure(
        list(
          sample_component = c(
            "0a45497c6bfbfb22_1", "0a45497c6bfbfb22_2",
            "0a45497c6bfbfb22_3", "0a45497c6bfbfb22_4",
            "0a45497c6bfbfb22_5", "0a45497c6bfbfb22_6"
          ),
          sample_alias = c(
            "S1", "S1", "S1", "S1",
            "S1", "S1"
          ),
          mean_coreness = c(
            2.5690926210872, 2.5690926210872,
            2.5690926210872, 2.5690926210872,
            2.5690926210872, 2.5690926210872
          ),
          percent_dangling_nodes = c(
            23.5629148198333, 23.5629148198333,
            23.5629148198333, 23.5629148198333,
            23.5629148198333, 23.5629148198333
          ),
          percent_well_connected_nodes = c(
            62.5519601313644, 62.5519601313644,
            62.5519601313644, 62.5519601313644,
            62.5519601313644, 62.5519601313644
          )
        ),
        class = c("grouped_df", "tbl_df", "tbl", "data.frame"),
        row.names = c(NA, -6L),
        groups = structure(
          list(
            sample_component = c(
              "0a45497c6bfbfb22_1", "0a45497c6bfbfb22_2",
              "0a45497c6bfbfb22_3", "0a45497c6bfbfb22_4",
              "0a45497c6bfbfb22_5", "0a45497c6bfbfb22_6"
            ),
            .rows = structure(list(1L, 2L, 3L, 4L, 5L, 6L),
                              ptype = integer(0),
                              class = c("vctrs_list_of", "vctrs_vctr", "list")
            )
          ),
          class = c("tbl_df", "tbl", "data.frame"),
          row.names = c(NA, -6L),
          .drop = TRUE
        )
      )
    )

    expect_equal(
      res$sample_summary,
      structure(
        list(
          sample_alias = "S1",
          median_mean_coreness = 2.47440345751157,
          median_percent_dangling_nodes = 24.993831729583,
          median_percent_well_connected_nodes = 57.1067689310933
        ),
        class = c("tbl_df", "tbl", "data.frame"),
        row.names = c(NA, -1L)
      )
    )


    expect_no_error(
      sample_qc_tables <-
        get_qc_metrics(seur, sample_qc_metrics, sample_sheet)
    )
    expect_no_error(
      tabl <- key_metric_table(sample_qc_tables)
    )
    expect_s3_class(tabl$sample, "datatables")

    expect_no_error(
      tabl <- key_metric_table(sample_qc_tables, return_data = TRUE)
    )
    expect_type(tabl, "list")
    expect_s3_class(tabl$sample, "tbl_df")
    expect_equal(
      tabl$sample,
      switch(data_type,
             default =
               structure(
                 list(
                   `Sample ID` = "S1",
                   `Number of cells` = "30",
                   `Number of cells >10k nodes` = "30",
                   `Median isotype % counts` = "0.04",
                   `Median intracellular % counts` = "0",
                   `Median proteins per cell [k]` = "48.64",
                   `Median reads per cell [k]` = "289.86",
                   `Q30 [%]` = "95.77",
                   `Total reads [M]` = "400.75",
                   `Valid reads [M]` = "359.95",
                   `Graph Nodes [M]` = "1.47",
                   `Graph Edges [M]` = "3.17",
                   `Graph node saturation [%]` = "83.2",
                   `Graph edge saturation [%]` = "63.74",
                   `Valid reads saturation [%]` = "35.4",
                   `Valid reads fraction [%]` = "89.82",
                   `Graph reads fraction [%]` = "2.18",
                   `% Denoised UMIs` = "2.67",
                   `Total denoised UMIs [M]` = "1.62",
                   `Median mean coreness` = "2.47",
                   `Median % dangling nodes` = "24.99",
                   `Median % well connected nodes` = "57.11",
                   `Top 3 % counts` = "35.2",
                   `Top 5 % counts` = "45.23",
                   `Top 5 markers` = "CD44, B2M, CD59, CD45, CD43",
                   `% Crossing edges (Initial)` = "5.06",
                   `% Crossing edges (Refinement)` = "0.01"
                 ),
                 row.names = c(NA, -1L),
                 class = c("tbl_df", "tbl", "data.frame")
               ),
             hashing =
               structure(list(
                 `Sample ID` = "S1", `Number of cells` = "30",
                 `Number of cells >10k nodes` = "30", `Median isotype % counts` = "0.04",
                 `Median intracellular % counts` = "0", `Median proteins per cell [k]` = "48.64",
                 `Median reads per cell [k]` = "289.86", `% Denoised UMIs` = "4.78",
                 `Total denoised UMIs [M]` = "0", `Median mean coreness` = "2.47",
                 `Median % dangling nodes` = "24.99", `Median % well connected nodes` = "57.11",
                 `Top 3 % counts` = "35.2", `Top 5 % counts` = "45.23", `Top 5 markers` = "CD44, B2M, CD59, CD45, CD43"
               ), row.names = c(
                 NA,
                 -1L
               ), class = c("tbl_df", "tbl", "data.frame"))
      )
    )
  })
}
