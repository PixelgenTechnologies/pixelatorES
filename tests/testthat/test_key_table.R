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
      switch(data_type,
        default =
          structure(list(sample_alias = c("S1", "S2"), top3_fraction = c(
            0.216648620555243,
            0.216648620555243
          ), top5_fraction = c(
            0.320659188358093,
            0.320659188358093
          ), top_markers = c(
            "B2M, HLA-DR-DP-DQ, CD45, CD45RA, HLA-ABC",
            "B2M, HLA-DR-DP-DQ, CD45, CD45RA, HLA-ABC"
          )), row.names = c(
            NA,
            -2L
          ), class = c("tbl_df", "tbl", "data.frame")),
        hashing =
          structure(list(sample_alias = c("S1", "S2", "S11", "S12"), top3_fraction = c(
            0.420525431861804,
            0.388412618367233, 0.39153322126869, 0.30195929854947
          ), top5_fraction = c(
            0.57017754318618,
            0.496273883565865, 0.539237370609573, 0.450205672223425
          ), top_markers = c(
            "B2M, CD45RA, CD45, HLA-ABC, CD43",
            "HLA-DR-DP-DQ, B2M, CD45RA, CD40, IgM", "B2M, CD45RA, CD45, HLA-ABC, CD43",
            "B2M, HLA-DR-DP-DQ, CD59, HLA-ABC, CD29"
          )), row.names = c(
            NA,
            -4L
          ), class = c("tbl_df", "tbl", "data.frame"))
      )
    )

    # Sequencing saturation
    expect_no_error(res <- get_seq_saturation(seur, sample_qc_metrics))
    expect_equal(
      res,
      switch(data_type,
        default =
          structure(list(
            sample_alias = c("S1", "S2"), Q30 = c(
              0.957683900982373,
              0.958092725597593
            ), total_reads = c(400745924L, 210603520L),
            deduped_valid_reads = c(
              232530522L,
              124914607L
            ), valid_reads = c(
              359952911L,
              188742157L
            ), graph_edges = c(
              384168L,
              384168L
            ), graph_proteins = c(
              230890L,
              230890L
            ), graph_reads = c(923084L, 923084L), valid_reads_saturation = c(
              35.3997384396761,
              33.8173257180694
            ), graph_edge_saturation = c(
              58.3821190704205,
              58.3821190704205
            ), graph_node_saturation = c(
              74.9871084321687,
              74.9871084321687
            ), fraction_valid_reads = c(
              89.8207291560625,
              89.6196592535585
            ), fraction_graph_reads = c(
              0.230341456947669,
              0.438304165096576
            )
          ), row.names = c(NA, -2L), class = c(
            "tbl_df",
            "tbl", "data.frame"
          )),
        hashing =
          structure(list(pool = c("pool1", "pool2"), total_reads = c(
            2240712L,
            2240712L
          ), deduped_valid_reads = 1159113:1159114, valid_reads = c(
            1953716L,
            1953716L
          ), graph_edges = c(66575L, 84144L), graph_proteins = c(
            40627L,
            48052L
          ), graph_reads = c(154863L, 203615L), valid_reads_saturation = c(
            40.6713667697864,
            40.6713155852744
          ), graph_edge_saturation = c(
            57.010389828429,
            58.6749502738011
          ), graph_node_saturation = c(
            73.7658446497872,
            76.400559880166
          ), fraction_valid_reads = c(
            87.1917497652532,
            87.1917497652532
          ), fraction_graph_reads = c(
            6.91132997011664,
            9.08706696799946
          )), row.names = c(NA, -2L), class = c(
            "tbl_df",
            "tbl", "data.frame"
          ))
      )
    )

    # Read stats
    expect_no_error(res <- get_read_stats(seur))
    expect_equal(
      res,
      switch(data_type,
        default =
          structure(list(sample_alias = structure(1:2, levels = c(
            "S1",
            "S2"
          ), class = "factor"), n_cells = c(24L, 24L), n_cells_over10k = c(
            7L,
            7L
          ), median_reads_per_cell = c(34696, 34696), median_abs_per_cell = c(
            8992.5,
            8992.5
          ), median_isotype_count_pct = c(
            0.190114870060899,
            0.190114870060899
          ), median_intracellular_count_pct = c(0, 0)), row.names = c(
            NA,
            -2L
          ), class = c("tbl_df", "tbl", "data.frame")),
        hashing =
          structure(list(sample_alias = structure(1:4, levels = c(
            "S1",
            "S2", "S11", "S12"
          ), class = "factor"), n_cells = c(
            2L, 2L, 2L,
            2L
          ), n_cells_over10k = c(1L, 1L, 1L, 1L), median_reads_per_cell = c(
            39277.5,
            38154, 56518.5, 45289
          ), median_abs_per_cell = c(
            9539.5, 10774,
            13326.5, 10699.5
          ), median_isotype_count_pct = c(
            0.1997935831379,
            0.299137580565269, 0.288901593366738, 0.197655770206966
          ), median_intracellular_count_pct = c(
            0,
            0, 0, 0
          )), row.names = c(NA, -4L), class = c(
            "tbl_df", "tbl",
            "data.frame"
          ))
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
      switch(data_type,
        default =
          structure(list(sample_alias = c(
            "S1", "S1", "S1", "S1", "S1",
            "S1"
          ), umi_type = c("umi1", "umi1", "umi1", "umi1", "umi1", "umi1"), degree = 1:6, n = c(
            27439568L, 7545696L, 5571202L, 4951939L,
            4609669L, 4290856L
          ), percent_nodes = c(
            19.0277655876346, 5.23250711102857,
            3.86330884281272, 3.43388549325068, 3.19654089191878, 2.97546237383531
          )), row.names = c(NA, -6L), class = c("tbl_df", "tbl", "data.frame")),
        hashing =
          structure(list(pool = c(
            "pool1", "pool1", "pool1", "pool1", "pool1",
            "pool1"
          ), umi_type = c(
            "umi1", "umi1", "umi1", "umi1", "umi1",
            "umi1"
          ), degree = 1:6, n = c(
            75201L, 55220L, 45622L, 33576L,
            22807L, 14510L
          ), percent_nodes = c(
            14.6206175184554, 10.7359011099467,
            8.86985295976078, 6.52786337681224, 4.43414879780071, 2.82104174402983
          )), row.names = c(NA, -6L), class = c("tbl_df", "tbl", "data.frame"))
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
      switch(data_type,
        default =
          structure(list(sample_component = c(
            "S1_128b348aca07cb57", "S1_128b348aca07cb57",
            "S1_128b348aca07cb57", "S1_128b348aca07cb57", "S1_2ee5115016f6d3bf",
            "S1_2ee5115016f6d3bf"
          ), sample_alias = structure(c(
            1L, 1L, 1L,
            1L, 1L, 1L
          ), levels = c("S1", "S2"), class = "factor"), molecules = c(
            8372L,
            8372L, 8372L, 8372L, 8457L, 8457L
          ), k_core = c(
            "k_core_1", "k_core_2",
            "k_core_3", "k_core_4", "k_core_1", "k_core_2"
          ), nodes = c(
            2033,
            3326, 3013, 0, 2312, 4187
          ), percent_nodes = c(
            24.2833253702819,
            39.7276636407071, 35.989010989011, 0, 27.3382996334398, 49.5092822513894
          ), coreness = c(1, 2, 3, 4, 1, 2)), row.names = c(NA, -6L), class = c(
            "tbl_df",
            "tbl", "data.frame"
          )),
        hashing =
          structure(list(
            sample_component = c(
              "S1_d2146defe08567d3", "S1_d2146defe08567d3",
              "S1_d2146defe08567d3", "S1_68189b2c75de4098", "S1_68189b2c75de4098",
              "S1_68189b2c75de4098"
            ), sample_alias = structure(c(
              1L, 1L, 1L,
              1L, 1L, 1L
            ), levels = c("S1", "S2", "S11", "S12"), class = "factor"),
            molecules = c(
              10562L, 10562L, 10562L, 8517L, 8517L, 8517L
            ), k_core = c(
              "k_core_1", "k_core_2", "k_core_3", "k_core_1",
              "k_core_2", "k_core_3"
            ), nodes = c(
              2693, 3849, 4020, 2150,
              3442, 2925
            ), percent_nodes = c(
              25.4970649498201, 36.4419617496686,
              38.0609733005113, 25.2436303862863, 40.413291064929, 34.3430785487848
            ), coreness = c(1, 2, 3, 1, 2, 3)
          ), row.names = c(NA, -6L), class = c("tbl_df", "tbl", "data.frame"))
      )
    )

    expect_equal(
      head(res$component_summary),
      switch(data_type,
        default =
          structure(list(sample_component = c(
            "S1_128b348aca07cb57", "S1_218e57362e5c1a3d",
            "S1_2ee5115016f6d3bf", "S1_39576a361fdfa5ca", "S1_500dadc305f7fd2d",
            "S1_6d7587035b64afd5"
          ), sample_alias = structure(c(
            1L, 1L, 1L,
            1L, 1L, 1L
          ), levels = c("S1", "S2"), class = "factor"), mean_coreness = c(
            2.11705685618729,
            1.98440613358746, 1.95814118481731, 2.08445612882024, 1.6720041322314,
            1.69893867924528
          ), percent_dangling_nodes = c(
            24.2833253702819,
            27.9130208784545, 27.3382996334398, 24.663161353927, 34.0650826446281,
            33.3254716981132
          ), percent_well_connected_nodes = c(
            35.989010989011,
            26.3536342372, 23.1524181151709, 33.1087742359514, 1.2654958677686,
            3.21933962264151
          )), row.names = c(NA, -6L), groups = structure(list(
            sample_component = c(
              "S1_128b348aca07cb57", "S1_218e57362e5c1a3d",
              "S1_2ee5115016f6d3bf", "S1_39576a361fdfa5ca", "S1_500dadc305f7fd2d",
              "S1_6d7587035b64afd5"
            ), .rows = structure(list(
              1L, 2L, 3L,
              4L, 5L, 6L
            ), ptype = integer(0), class = c(
              "vctrs_list_of",
              "vctrs_vctr", "list"
            ))
          ), row.names = c(NA, -6L), class = c(
            "tbl_df",
            "tbl", "data.frame"
          ), .drop = TRUE), class = c(
            "grouped_df",
            "tbl_df", "tbl", "data.frame"
          )),
        hashing =
          structure(list(
            sample_component = c(
              "S11_19b04397ed7f04ba", "S11_fe556695f452a4bb",
              "S12_3d23c6539cbead8d", "S12_63fab986007319b4", "S1_68189b2c75de4098",
              "S1_d2146defe08567d3"
            ), sample_alias = structure(c(
              3L, 3L, 4L,
              4L, 1L, 1L
            ), levels = c("S1", "S2", "S11", "S12"), class = "factor"),
            mean_coreness = c(
              2.19097736117963, 2.08216409036861, 2.30304600082884,
              1.77577253766919, 2.09099448162499, 2.12563908350691
            ), percent_dangling_nodes = c(
              23.0170476347092,
              24.1617122473246, 20.6589307915458, 30.4843789903805, 25.2436303862863,
              25.4970649498201
            ), percent_well_connected_nodes = c(
              42.1147837526723,
              32.3781212841855, 50.9635308744302, 8.06163275729974, 34.3430785487848,
              38.0609733005113
            )
          ), row.names = c(NA, -6L), groups = structure(list(
            sample_component = c(
              "S11_19b04397ed7f04ba", "S11_fe556695f452a4bb",
              "S12_3d23c6539cbead8d", "S12_63fab986007319b4", "S1_68189b2c75de4098",
              "S1_d2146defe08567d3"
            ), .rows = structure(list(
              1L, 2L, 3L,
              4L, 5L, 6L
            ), ptype = integer(0), class = c(
              "vctrs_list_of",
              "vctrs_vctr", "list"
            ))
          ), row.names = c(NA, -6L), class = c(
            "tbl_df",
            "tbl", "data.frame"
          ), .drop = TRUE), class = c(
            "grouped_df",
            "tbl_df", "tbl", "data.frame"
          ))
      )
    )

    expect_equal(
      res$sample_summary,
      switch(data_type,
        default =
          structure(list(sample_alias = structure(1:2, levels = c(
            "S1",
            "S2"
          ), class = "factor"), median_mean_coreness = c(
            1.95019805561644,
            1.95019805561644
          ), median_percent_dangling_nodes = c(
            26.1241628830854,
            26.1241628830854
          ), median_percent_well_connected_nodes = c(
            22.804941443752,
            22.804941443752
          )), row.names = c(NA, -2L), class = c(
            "tbl_df",
            "tbl", "data.frame"
          )),
        hashing =
          structure(list(sample_alias = structure(1:4, levels = c(
            "S1",
            "S2", "S11", "S12"
          ), class = "factor"), median_mean_coreness = c(
            2.10831678256595,
            1.83892332963412, 2.13657072577412, 2.03940926924902
          ), median_percent_dangling_nodes = c(
            25.3703476680532,
            29.6128537644793, 23.5893799410169, 25.5716548909632
          ), median_percent_well_connected_nodes = c(
            36.202025924648,
            13.5051867278914, 37.2464525184289, 29.512581815865
          )), row.names = c(
            NA,
            -4L
          ), class = c("tbl_df", "tbl", "data.frame"))
      )
    )

    # Hash stats
    if (data_type == "hashing") {
      expect_no_error(res <- get_hash_stats(seur, sample_qc_metrics))

      expect_equal(
        lapply(res[-6], head),
        list(pool_stats = structure(list(
          pool = c("pool1", "pool2"),
          cells_in_pool = c(216L, 290L), percent_called_cells = c(
            81.4814814814815,
            82.7586206896552
          )
        ), row.names = c(NA, -2L), class = c(
          "tbl_df",
          "tbl", "data.frame"
        )), component_stats = structure(list(
          component = c(
            "S1_d2146defe08567d3",
            "S1_d2146defe08567d3", "S1_68189b2c75de4098", "S1_68189b2c75de4098",
            "S11_19b04397ed7f04ba", "S11_19b04397ed7f04ba"
          ), sample_alias = structure(c(
            1L,
            1L, 1L, 1L, 3L, 3L
          ), levels = c("S1", "S2", "S11", "S12"), class = "factor"),
          id = c("CD298", "CD298", "CD298", "CD298", "CD298", "CD298"), version = c("1", "2", "1", "2", "1", "2"), count = c(
            0,
            0, 1, 0, 0, 0
          ), purity = c(0, 0, 1, 0, 0, 0)
        ), row.names = c(
          NA,
          -6L
        ), groups = structure(list(
          component = c(
            "S11_19b04397ed7f04ba",
            "S1_68189b2c75de4098", "S1_d2146defe08567d3"
          ), sample_alias = structure(c(
            3L,
            1L, 1L
          ), levels = c("S1", "S2", "S11", "S12"), class = "factor"),
          id = c("CD298", "CD298", "CD298"), .rows = structure(list(
            5:6, 3:4, 1:2
          ), ptype = integer(0), class = c(
            "vctrs_list_of",
            "vctrs_vctr", "list"
          ))
        ), row.names = c(NA, -3L), class = c(
          "tbl_df",
          "tbl", "data.frame"
        ), .drop = TRUE), class = c(
          "grouped_df",
          "tbl_df", "tbl", "data.frame"
        )), component_stats_heatmap_purity = structure(list(
          sample_component = c(
            "S1_d2146defe08567d3", "S1_68189b2c75de4098",
            "S11_19b04397ed7f04ba", "S11_fe556695f452a4bb", "S12_3d23c6539cbead8d",
            "S12_63fab986007319b4"
          ), sample_alias = structure(c(
            1L, 1L,
            3L, 3L, 4L, 4L
          ), levels = c("S1", "S2", "S11", "S12"), class = "factor"),
          `CD298-1` = c(0, 1, 0, 1, 0, 0), `CD298-2` = c(
            0, 0, 0, 0,
            1, 1
          ), `B2M-1` = c(
            0.964071856287425, 0.952380952380952,
            0.967153284671533, 0.990950226244344, 0, 0
          ), `B2M-2` = c(
            0.0359281437125748,
            0.0476190476190476, 0.0328467153284672, 0.00904977375565611,
            1, 1
          ), `CD98-1` = c(
            0.666666666666667, 1, 0.666666666666667,
            1, 0, 0
          ), `CD98-2` = c(
            0.333333333333333, 0, 0.333333333333333,
            0, 0, 1
          ), molecules = c(
            10562L, 8517L, 18243L, 8410L, 9652L,
            11747L
          )
        ), row.names = c(NA, -6L), class = c(
          "tbl_df", "tbl",
          "data.frame"
        ), id_hash_order = c(
          "CD298-1", "B2M-1", "CD98-1",
          "CD298-2", "B2M-2", "CD98-2"
        )), component_stats_heatmap_fraction = structure(list(
          sample_component = c(
            "S1_d2146defe08567d3", "S1_68189b2c75de4098",
            "S11_19b04397ed7f04ba", "S11_fe556695f452a4bb", "S12_3d23c6539cbead8d",
            "S12_63fab986007319b4"
          ), sample_alias = structure(c(
            1L, 1L,
            3L, 3L, 4L, 4L
          ), levels = c("S1", "S2", "S11", "S12"), class = "factor"),
          `CD298-1` = c(
            0, 0.0112359550561798, 0, 0.0131578947368421,
            0, 0
          ), `CD298-2` = c(0, 0, 0, 0, 0.4, 0.0611764705882353),
          `B2M-1` = c(
            0.930635838150289, 0.898876404494382, 0.946428571428571,
            0.960526315789474, 0, 0
          ), `B2M-2` = c(
            0.0346820809248555,
            0.0449438202247191, 0.0321428571428571, 0.0087719298245614,
            0.6, 0.68
          ), `CD98-1` = c(
            0.023121387283237, 0.0449438202247191,
            0.0142857142857143, 0.0175438596491228, 0, 0
          ), `CD98-2` = c(
            0.0115606936416185,
            0, 0.00714285714285714, 0, 0, 0.258823529411765
          )
        ), row.names = c(
          NA,
          -6L
        ), groups = structure(list(.rows = structure(list(
          1L, 2L,
          3L, 4L, 5L, 6L
        ), ptype = integer(0), class = c(
          "vctrs_list_of",
          "vctrs_vctr", "list"
        ))), row.names = c(NA, -6L), class = c(
          "tbl_df",
          "tbl", "data.frame"
        )), id_hash_order = c(
          "CD298-1", "B2M-1",
          "CD98-1", "CD298-2", "B2M-2", "CD98-2"
        ), class = c(
          "rowwise_df",
          "tbl_df", "tbl", "data.frame"
        )), sample_stats = structure(list(
          sample_alias = structure(c(1L, 4L, 2L, 3L), levels = c(
            "S1",
            "S11", "S12", "S2"
          ), class = "factor"), mean_purity_B2M = c(
            0.958226404334189,
            0.996124031007752, 0.979051755457938, 1
          ), mean_purity_CD298 = c(
            0.5,
            1, 0.5, 1
          ), mean_purity_CD98 = c(
            0.833333333333333, 1, 0.833333333333333,
            0.5
          ), hash_pct = c(
            0.0137323759106871, 0.0293762762205309,
            0.0190597681311672, 0.0200943969344362
          )
        ), row.names = c(
          NA,
          -4L
        ), class = c("tbl_df", "tbl", "data.frame")))
      )

      expect_equal(
        head(res[[6]]),
        structure(list(pool = c(
          "pool1", "pool1", "pool1", "pool1", "pool1",
          "pool1"
        ), sample_alias = c("S1", "S1", "S1", "S1", "S1", "S1"), sample_confidence = c(
          0.988843436221644, 0.987964696442899,
          0.96164199192463, 0.994916690200508, 0.982078853046595, 0.988266464799394
        )), row.names = c(NA, -6L), class = c("tbl_df", "tbl", "data.frame"))
      )
    }


    # Key table
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
          structure(list(
            `Sample ID` = c("S1", "S2"), `Number of cells in sample` = c(
              "24",
              "24"
            ), `Number of cells >10k nodes` = c("7", "7"), `Median isotype % counts` = c(
              "0.19",
              "0.19"
            ), `Median intracellular % counts` = c("0", "0"), `Median proteins per cell [k]` = c(
              "8.99",
              "8.99"
            ), `Median reads per cell [k]` = c("34.7", "34.7"), `Q30 [%]` = c(
              "95.77",
              "95.81"
            ), `Total reads [M]` = c("400.75", "210.6"), `Valid reads [M]` = c(
              "359.95",
              "188.74"
            ), `Graph Nodes [M]` = c("0.23", "0.23"), `Graph Edges [M]` = c(
              "0.38",
              "0.38"
            ), `Graph node saturation [%]` = c("74.99", "74.99"), `Graph edge saturation [%]` = c(
              "58.38",
              "58.38"
            ), `Valid reads saturation [%]` = c("35.4", "33.82"),
            `Valid reads fraction [%]` = c("89.82", "89.62"), `Graph reads fraction [%]` = c(
              "0.23",
              "0.44"
            ), `% Denoised UMIs` = c("2.67", "4.86"), `Total denoised UMIs [M]` = c(
              "1.62",
              "2.72"
            ), `Median mean coreness` = c("1.95", "1.95"), `Median % dangling nodes` = c(
              "26.12",
              "26.12"
            ), `Median % well connected nodes` = c("22.8", "22.8"),
            `Top 3 % counts` = c("21.66", "21.66"), `Top 5 % counts` = c(
              "32.07",
              "32.07"
            ), `Top 5 markers` = c(
              "B2M, HLA-DR-DP-DQ, CD45, CD45RA, HLA-ABC",
              "B2M, HLA-DR-DP-DQ, CD45, CD45RA, HLA-ABC"
            ), `% Crossing edges (Initial)` = c(
              "5.06",
              "4.27"
            ), `% Crossing edges (Refinement)` = c("0.01", "0.02")
          ), row.names = c(NA, -2L), class = c("tbl_df", "tbl", "data.frame")),
        hashing =
          structure(list(
            `Sample ID` = c("S1", "S2", "S11", "S12"), `Number of cells in sample` = c(
              "2",
              "2", "2", "2"
            ), `Number of cells >10k nodes` = c(
              "1", "1", "1",
              "1"
            ), `Median isotype % counts` = c("0.2", "0.3", "0.29", "0.2"),
            `Median intracellular % counts` = c("0", "0", "0", "0"), `Median proteins per cell [k]` = c(
              "9.54",
              "10.77", "13.33", "10.7"
            ), `Median reads per cell [k]` = c(
              "39.28",
              "38.15", "56.52", "45.29"
            ), `% Denoised UMIs` = c(
              "4.78", "4.61",
              "4.79", "4.73"
            ), `Total denoised UMIs [M]` = c(
              "0", "0.01", "0",
              "0.01"
            ), `Median mean coreness` = c("2.11", "1.84", "2.14", "2.04"),
            `Median % dangling nodes` = c(
              "25.37", "29.61", "23.59", "25.57"
            ), `Median % well connected nodes` = c(
              "36.2", "13.51", "37.25",
              "29.51"
            ), `Top 3 % counts` = c("42.05", "38.84", "39.15", "30.2"),
            `Top 5 % counts` = c("57.02", "49.63", "53.92", "45.02"),
            `Top 5 markers` = c(
              "B2M, CD45RA, CD45, HLA-ABC, CD43", "HLA-DR-DP-DQ, B2M, CD45RA, CD40, IgM",
              "B2M, CD45RA, CD45, HLA-ABC, CD43", "B2M, HLA-DR-DP-DQ, CD59, HLA-ABC, CD29"
            ), `% B2M hash purity` = c("95.82", "99.61", "97.91", "100"),
            `% CD298 hash purity` = c("50", "100", "50", "100"), `% CD98 hash purity` = c(
              "83.33",
              "100", "83.33", "50"
            ), `% hash counts` = c(
              "1.37", "2.94",
              "1.91", "2.01"
            )
          ), row.names = c(NA, -4L), class = c(
            "tbl_df",
            "tbl", "data.frame"
          ))
      )
    )
  })
}
