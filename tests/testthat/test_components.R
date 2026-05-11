library(dplyr)

data_types <- c("default", "hashing")

for (data_type in data_types) {
  pg_data <- get_test_data(type = data_type)
  pg_data_small <- get_test_data(type = data_type)
  sample_qc_metrics <- get_test_qc_metrics(type = data_type)
  sample_sheet <- read_samplesheet(test_samplesheet(type = data_type))
  qc_metrics_tables <-
    get_qc_metrics(pg_data, sample_qc_metrics, sample_sheet)


  test_message <- paste("Components work as expected for", data_type, "data")

  test_that(test_message, {
    # component_control_markers
    expect_no_error(
      component <- component_control_markers(pg_data)
    )

    expect_s3_class(component$p1, "ggplot")
    expect_no_error(
      ggplot2::ggplot_build(component$p1)
    )
    expect_s3_class(component$p2, "ggplot")
    expect_no_error(
      ggplot2::ggplot_build(component$p2)
    )
    expect_s3_class(component$tabl, "datatables")

    # component_crossing_edges
    expect_no_error(component <- component_crossing_edges(qc_metrics_tables))
    expect_no_error(
      ggplot2::ggplot_build(component$plot)
    )

    # component_cell_recovery
    expect_no_error(component <- component_cell_recovery(pg_data, sample_qc_metrics, sample_levels = NULL))
    expect_no_error(
      ggplot2::ggplot_build(component$plots[[1]])
    )

    for (plot in component$plots[[2]]) {
      if (data_type == "hashing") {
        for (subplot in plot$Pools) {
          expect_s3_class(subplot, "ggplot")
          expect_no_error(
            ggplot2::ggplot_build(subplot)
          )
        }

        for (subplot in plot$Samples) {
          expect_s3_class(subplot, "ggplot")
          expect_no_error(
            ggplot2::ggplot_build(subplot)
          )
        }
      } else {
        expect_s3_class(plot, "ggplot")
        expect_no_error(
          ggplot2::ggplot_build(plot)
        )
      }
    }

    for (table in component$table) {
      expect_s3_class(table, "datatables")
    }

    # component_node_degree
    expect_no_error(
      component <- component_node_degree(pg_data, sample_levels = NULL)
    )

    expect_s3_class(component$plot, "ggplot")
    expect_no_error(
      ggplot2::ggplot_build(component$plot)
    )
    expect_s3_class(component$table, "datatables")

    # component_node_edge_count
    expect_no_error(
      component <- component_node_edge_count(sample_qc_metrics, sample_levels = NULL)
    )

    for (plot in component$plots) {
      expect_s3_class(plot, "ggplot")
      expect_no_error(
        ggplot2::ggplot_build(plot)
      )
    }
    expect_s3_class(component$table, "datatables")

    # component_sequencing_saturation
    expect_no_error(
      component <- component_sequencing_saturation(qc_metrics_tables, sample_levels = NULL)
    )

    for (plot in component$plots) {
      expect_s3_class(plot, "ggplot")
      expect_no_error(
        ggplot2::ggplot_build(plot)
      )
    }
    expect_s3_class(component$table, "datatables")

    # component_sequencing_reads_per_cell
    expect_no_error(component <- component_sequencing_reads_per_cell(pg_data))
    expect_s3_class(component$plot, "ggplot")
    expect_no_error(
      ggplot2::ggplot_build(component$plot)
    )
    expect_s3_class(component$table, "datatables")

    # component_sequencing_reads_and_molecules
    expect_no_error(component <- component_sequencing_reads_and_molecules(sample_qc_metrics, sample_levels = NULL))

    for (plot in component$plots) {
      expect_s3_class(plot, "ggplot")
      expect_no_error(
        ggplot2::ggplot_build(plot)
      )
    }
    expect_s3_class(component$table, "datatables")

    # component_abundance_per_celltype
    temp <-
      pg_data %>%
      subset(features = rownames(pg_data)[1:3])

    temp[["celltype"]] <- sample(
      c("CD4 T", "CD8 T", "B"),
      ncol(temp),
      replace = TRUE
    )
    temp[["seurat_clusters"]] <- 1

    expect_no_error(
      component <- component_abundance_per_celltype(
        temp,
        params = list(
          control_markers = c("mIgG1", "mIgG2a", "mIgG2b")
        ),
        sample_palette = c("red", "black", "blue")
      )
    )

    expect_s3_class(component[[1]], "ggplot")
    expect_no_error(
      ggplot2::ggplot_build(component[[1]])
    )
    # component_proximity_per_marker
    temp <-
      pg_data_small %>%
      subset(features = rownames(pg_data_small)[1:3])

    set.seed(37)
    temp[["celltype"]] <- sample(
      c("T", "B"),
      ncol(temp),
      replace = TRUE
    )
    temp[["seurat_clusters"]] <- 1
    temp[["condition"]] <- "good"

    proximity_scores <-
      filter_proximity_scores(
        temp,
        list(
          control_markers = c("mIgG1", "mIgG2a", "mIgG2b")
        )
      )

    expect_no_error(
      component <- component_proximity_per_marker(
        proximity_scores,
        sample_palette = c("red", "black")
      )
    )
    expect_s3_class(component[[1]], "ggplot")
    expect_no_error(
      ggplot2::ggplot_build(component[[1]])
    )
    expect_equal(
      component[[1]]$data,
      switch(data_type,
        default =
          structure(list(
            sample_alias = structure(c(
              1L, 1L, 1L, 1L, 1L,
              1L, 1L, 2L, 2L, 2L, 2L, 2L, 2L, 2L
            ), levels = c("S1", "S2"), class = "factor"),
            celltype = c(
              "B", "B", "Mono & DC", "NK", "Platelets", "T",
              "T", "B", "B", "B", "Mono & DC", "NK", "Platelets", "T"
            ),
            marker_1 = structure(c(
              1L, 1L, NA, NA, NA, 1L, 1L, 1L, 1L,
              1L, NA, NA, NA, 1L
            ), levels = c("CD11b", "B2M", "HLA-ABC"), class = "factor"), marker_2 = structure(c(
              1L, 1L, NA,
              NA, NA, 1L, 1L, 1L, 1L, 1L, NA, NA, NA, 1L
            ), levels = c(
              "CD11b",
              "B2M", "HLA-ABC"
            ), class = "factor"), join_count = c(
              0, 0,
              NA, NA, NA, 0, 0, 0, 0, 0, NA, NA, NA, 0
            ), join_count_expected_mean = c(
              0.2,
              0.12, NA, NA, NA, 0.04, 0, 0.2, 0, 0.12, NA, NA, NA, 0.04
            ), join_count_expected_sd = c(
              0.449466574975495, 0.32659863237109,
              NA, NA, NA, 0.196946385566932, 0, 0.449466574975495, 0, 0.32659863237109,
              NA, NA, NA, 0.196946385566932
            ), join_count_z = c(
              -0.2, -0.12,
              NA, NA, NA, -0.04, 0, -0.2, 0, -0.12, NA, NA, NA, -0.04
            ),
            join_count_p = c(
              0.420740290560897, 0.452241573979416, NA,
              NA, NA, 0.484046563147169, 0.5, 0.420740290560897, 0.5, 0.452241573979416,
              NA, NA, NA, 0.484046563147169
            ), log2_ratio = c(
              0, 0, NA,
              NA, NA, 0, 0, 0, 0, 0, NA, NA, NA, 0
            ), sample_component = c(
              "S1_e2055911bf1693bd",
              "S1_dc8d14e01732603e", NA, NA, NA, "S1_99aa508551f33cd3",
              "S1_2ee5115016f6d3bf", "S2_e2055911bf1693bd", "S2_2ee5115016f6d3bf",
              "S2_dc8d14e01732603e", NA, NA, NA, "S2_99aa508551f33cd3"
            ),
            count_1 = c(
              49L, 21L, NA, NA, NA, 19L, 12L, 49L, 12L, 21L,
              NA, NA, NA, 19L
            ), count_2 = c(
              49L, 21L, NA, NA, NA, 19L,
              12L, 49L, 12L, 21L, NA, NA, NA, 19L
            ), p1 = c(
              0.016327890703099,
              0.010989010989011, NA, NA, NA, 0.00842945874001775, 0.0184331797235023,
              0.016327890703099, 0.0184331797235023, 0.010989010989011,
              NA, NA, NA, 0.00842945874001775
            ), p2 = c(
              0.016327890703099,
              0.010989010989011, NA, NA, NA, 0.00842945874001775, 0.0184331797235023,
              0.016327890703099, 0.0184331797235023, 0.010989010989011,
              NA, NA, NA, 0.00842945874001775
            ), condition = c(
              "good", "good",
              NA, NA, NA, "good", "good", "good", "good", "good", NA, NA,
              NA, "good"
            ), seurat_clusters = c(
              "1", "1", NA, NA, NA, "1",
              "1", "1", "1", "1", NA, NA, NA, "1"
            )
          ), row.names = c(
            NA,
            -14L
          ), class = c("tbl_df", "tbl", "data.frame")),
        hashing =
          structure(list(
            sample_alias = structure(c(
              1L, 1L, 1L, 1L, 1L,
              2L, 2L, 2L, 2L, 2L, 2L, 3L, 3L, 3L, 3L, 3L, 4L, 4L, 4L, 4L, 4L
            ), levels = c("S1", "S11", "S12", "S2"), class = "factor"), celltype = c(
              "B",
              "Mono & DC", "NK", "Platelets", "T", "B", "B", "Mono & DC", "NK",
              "Platelets", "T", "B", "Mono & DC", "NK", "Platelets", "T", "B",
              "Mono & DC", "NK", "Platelets", "T"
            ), marker_1 = structure(c(
              1L,
              NA, NA, NA, 1L, 1L, 1L, NA, NA, NA, NA, 1L, NA, NA, NA, NA, 1L,
              NA, NA, NA, 1L
            ), levels = c("CD45", "B2M", "HLA-DR-DP-DQ"), class = "factor"),
            marker_2 = structure(c(
              1L, NA, NA, NA, 1L, 1L, 1L, NA, NA,
              NA, NA, 1L, NA, NA, NA, NA, 1L, NA, NA, NA, 1L
            ), levels = c(
              "CD45",
              "B2M", "HLA-DR-DP-DQ"
            ), class = "factor"), join_count = c(
              153,
              NA, NA, NA, 345, 214, 258, NA, NA, NA, NA, 17, NA, NA, NA,
              NA, 1, NA, NA, NA, 18
            ), join_count_expected_mean = c(
              143.82,
              NA, NA, NA, 325.28, 162.28, 245.94, NA, NA, NA, NA, 15.12,
              NA, NA, NA, NA, 1.08, NA, NA, NA, 19.49
            ), join_count_expected_sd = c(
              12.001161559944,
              NA, NA, NA, 19.0899412934105, 13.1978571875678, 16.6192173052278,
              NA, NA, NA, NA, 3.47365263465745, NA, NA, NA, NA, 0.849004170076968,
              NA, NA, NA, 4.54938113229357
            ), join_count_z = c(
              0.764925957720613,
              NA, NA, NA, 1.03300474825489, 3.91881797665758, 0.725665943137187,
              NA, NA, NA, NA, 0.541217040887393, NA, NA, NA, NA, -0.0800000000000001,
              NA, NA, NA, -0.327517074668311
            ), join_count_p = c(
              0.222157817923894,
              NA, NA, NA, 0.150800838294423, 4.44921405915763e-05, 0.234021792365565,
              NA, NA, NA, NA, 0.294178996582506, NA, NA, NA, NA, 0.468118627986013,
              NA, NA, NA, 0.371638415373185
            ), log2_ratio = c(
              0.0892673380970873,
              NA, NA, NA, 0.0849142415955236, 0.399125588969599, 0.0690646698420441,
              NA, NA, NA, NA, 0.169076606803992, NA, NA, NA, NA, -0.111031312388744,
              NA, NA, NA, -0.114737184040853
            ), sample_component = c(
              "S1_d2146defe08567d3",
              NA, NA, NA, "S1_68189b2c75de4098", "S11_fe556695f452a4bb",
              "S11_19b04397ed7f04ba", NA, NA, NA, NA, "S12_63fab986007319b4",
              NA, NA, NA, NA, "S2_5bfdf506169806f2", NA, NA, NA, "S2_85f5bcdc05a7b286"
            ), count_1 = c(
              934L, NA, NA, NA, 1274L, 906L, 1582L, NA,
              NA, NA, NA, 345L, NA, NA, NA, NA, 77L, NA, NA, NA, 409L
            ),
            count_2 = c(
              934L, NA, NA, NA, 1274L, 906L, 1582L, NA, NA,
              NA, NA, 345L, NA, NA, NA, NA, 77L, NA, NA, NA, 409L
            ), p1 = c(
              0.365271802894016,
              NA, NA, NA, 0.594493700419972, 0.46749226006192, 0.366373320981936,
              NA, NA, NA, NA, 0.109211775878443, NA, NA, NA, NA, 0.0275985663082437,
              NA, NA, NA, 0.116690442225392
            ), p2 = c(
              0.365271802894016,
              NA, NA, NA, 0.594493700419972, 0.46749226006192, 0.366373320981936,
              NA, NA, NA, NA, 0.109211775878443, NA, NA, NA, NA, 0.0275985663082437,
              NA, NA, NA, 0.116690442225392
            ), condition = c(
              "good", NA,
              NA, NA, "good", "good", "good", NA, NA, NA, NA, "good", NA,
              NA, NA, NA, "good", NA, NA, NA, "good"
            ), seurat_clusters = c(
              "1",
              NA, NA, NA, "1", "1", "1", NA, NA, NA, NA, "1", NA, NA, NA,
              NA, "1", NA, NA, NA, "1"
            )
          ), row.names = c(NA, -21L), class = c(
            "tbl_df",
            "tbl", "data.frame"
          ))
      )
    )

    # component_proximity_selected
    set.seed(37)
    expect_no_error(
      component <- component_proximity_selected(
        temp,
        proximity_scores %>%
          filter(as.character(marker_1) == as.character(marker_2)),
        sample_palette = c("red", "black"),
        selected_contrasts = FALSE,
        proximity_score = "log2_ratio"
      )
    )

    expect_named(component, expected = switch(data_type,
      default = c("B2M/B2M", "CD11b/CD11b", "HLA-ABC/HLA-ABC"),
      hashing = c("B2M/B2M", "CD45/CD45", "HLA-DR-DP-DQ/HLA-DR-DP-DQ")
    ))
    expect_s3_class(component[[1]], "ggplot")
    expect_no_error(
      ggplot2::ggplot_build(component[[1]])
    )


    # component_dimred_plots
    expect_no_error(
      component <- component_dimred_plots(pg_data, sample_palette = c("red", "black", "blue"))
    )

    expect_s3_class(component$PCA, "ggplot")

    # component_proximity_heatmap_sample
    expect_no_error(
      component <- component_proximity_heatmap_sample(
        proximity_scores,
        temp,
        heatmap_gradient = c("red", "cyan"),
        n_markers = 2,
        plot_markers = NULL,
        test_mode = TRUE
      )
    )
    # component_proximity_heatmap_celltype
    expect_no_error(
      component <- component_proximity_heatmap_celltype(
        proximity_scores,
        temp,
        heatmap_gradient = c("red", "cyan"),
        n_markers = 2,
        plot_markers = NULL,
        test_mode = TRUE
      )
    )
  })
}
