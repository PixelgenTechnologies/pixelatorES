library(dplyr)

data_types <- c("default", "hashing")

cherry_gradient <-
  PixelgenGradient(100, "Cherry")

for (data_type in data_types) {
  # Read in Seurat objects
  pg_data <- get_test_data(type = data_type)
  pg_data_small <- get_test_data(type = data_type)

  # Get data folder and samplesheet
  data_folder <- test_data_folder(type = data_type)
  sample_sheet <- read_samplesheet(test_samplesheet(type = data_type))

  # Get QC metrics
  file_paths <-
    get_file_paths(
      data_folder = data_folder,
      sample_sheet = sample_sheet,
      stages = get_es_workflow_stages("amplicon_demux")
    )

  sample_qc_metrics <- read_qc_files(file_paths, sample_sheet)


  qc_metrics_tables <-
    get_qc_metrics(pg_data, sample_qc_metrics, sample_sheet)

  es_data <- test_es_data(
    params = list(
      control_markers = c("mIgG1", "mIgG2a", "mIgG2b")
    ),
    samplesheet = sample_sheet,
    effective_samplesheet = sample_sheet,
    file_paths = file_paths,
    pxl_data_processed = pg_data,
    qc_raw = sample_qc_metrics,
    qc = qc_metrics_tables,
    proximity = NULL
  )

  test_message <- paste("Components work as expected for", data_type, "data")

  test_that(test_message, {
    # component_control_markers
    expect_no_error(
      component <- component_control_markers(es_data)
    )

    expect_s3_class(component$p1, "ggplot")
    expect_equal(component$p1$labels$title, "Isotype control marker fraction")
    expect_equal(component$p1$labels$y, "Isotype control marker fraction [%]")
    expect_no_error(
      ggplot2::ggplot_build(component$p1)
    )
    expect_s3_class(component$p2, "ggplot")
    expect_equal(component$p2$labels$title, "Isotype control marker counts")
    expect_equal(component$p2$labels$y, "Isotype control marker counts")
    expect_no_error(
      ggplot2::ggplot_build(component$p2)
    )
    expect_s3_class(component$tabl, "datatables")

    # component_crossing_edges
    expect_no_error(component <- component_crossing_edges(es_data))
    expect_no_error(
      ggplot2::ggplot_build(component$plot)
    )

    # component_cell_recovery
    expect_no_error(component <- component_cell_recovery(es_data))
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

    # Molecule rank plots group components by exclusion and sample in focus
    MRP_plot <-
      if (data_type == "hashing") {
        component$plots[[2]]$Pools[[1]]
      } else {
        component$plots[[2]][[1]]
      }

    expect_equal(
      levels(MRP_plot$data$point_group),
      c(
        "Excluded, other samples",
        "Excluded, selected sample",
        "Included, other samples",
        "Included, selected sample"
      )
    )
    expect_equal(
      MRP_plot$labels$subtitle,
      "Excluded components shown faded"
    )
    # One named plot per sample: the plot lists are built from the sample_alias
    # factor levels, so a demoted factor would silently empty them.
    MRP_plot_lists <-
      if (data_type == "hashing") {
        list(component$plots[[2]]$Pools, component$plots[[2]]$Samples)
      } else {
        list(component$plots[[2]])
      }

    for (plot_list in MRP_plot_lists) {
      expect_gt(length(plot_list), 0)
      expect_true(all(nzchar(names(plot_list))))
    }

    if (data_type == "hashing") {
      sample_plot <- component$plots[[2]]$Samples[[1]]
      expect_true(any(!is.na(sample_plot$data$min_size_theshold)))
    }

    for (table in component$table) {
      expect_s3_class(table, "datatables")
    }

    # component_node_degree
    expect_no_error(
      component <- component_node_degree(es_data)
    )

    expect_s3_class(component$plot, "ggplot")
    expect_no_error(
      ggplot2::ggplot_build(component$plot)
    )
    expect_s3_class(component$table, "datatables")

    # component_node_edge_count
    expect_no_error(
      component <- component_node_edge_count(es_data)
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
      component <- component_sequencing_saturation(es_data)
    )

    for (plot in component$plots) {
      expect_s3_class(plot, "ggplot")
      expect_no_error(
        ggplot2::ggplot_build(plot)
      )
    }
    expect_s3_class(component$table, "datatables")

    # component_sequencing_reads_per_cell
    expect_no_error(component <- component_sequencing_reads_per_cell(es_data))
    expect_s3_class(component$plot, "ggplot")
    expect_no_error(
      ggplot2::ggplot_build(component$plot)
    )
    expect_s3_class(component$table, "datatables")

    # component_sequencing_reads_and_molecules
    expect_no_error(component <- component_sequencing_reads_and_molecules(es_data))

    for (plot in component$plots) {
      expect_s3_class(plot, "ggplot")
      expect_no_error(
        ggplot2::ggplot_build(plot)
      )
    }
    expect_s3_class(component$table, "datatables")

    # component_sequencing_saturation_curve
    expect_no_error(
      component <- component_sequencing_saturation_curve(es_data)
    )

    # component_denoising
    expect_no_error(component <- component_denoising(es_data))
    expect_type(component, "list")
    expect_named(component, c("plots", "tables"))
    expect_named(component$plots, c("removed_umis", "by_method", "isotype_reduction"))
    expect_named(component$tables, c("removed_umis", "by_method", "isotype_reduction"))

    for (plot in component$plots) {
      expect_s3_class(plot, "ggplot")
      expect_no_error(ggplot2::ggplot_build(plot))
    }
    for (table in component$tables) {
      expect_s3_class(table, "datatables")
    }

    temp <- qc_metrics_tables
    temp$denoising <-
      temp$denoising %>%
      rename(pool = 1)

    temp_es_data <- es_data
    temp_es_data$qc <- temp
    expect_no_error(component <- component_denoising(temp_es_data))
    expect_s3_class(component$plots$removed_umis, "ggplot")
    expect_no_error(
      ggplot2::ggplot_build(component$plots$removed_umis)
    )

    # component_abundance_per_marker
    temp <-
      pg_data %>%
      subset(features = rownames(pg_data)[1:3])

    temp[["celltype"]] <- sample(
      c("CD4 T", "CD8 T", "B"),
      ncol(temp),
      replace = TRUE
    )
    temp[["seurat_clusters"]] <- 1

    abundance_es_data <- es_data
    abundance_es_data$pxl_data_processed <- temp
    expect_no_error(
      component <- component_abundance_per_marker(
        abundance_es_data,
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

    proximity_es_data <- es_data
    proximity_es_data$proximity <- proximity_scores
    expect_no_error(
      component <- component_proximity_per_marker(
        proximity_es_data,
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
          structure(list(sample_alias = structure(c(1L, 1L, 1L, 1L, 1L,
          1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L,
          1L, 1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L,
          2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L), levels = c("S1",
          "S2"), class = "factor"), celltype = c("B", "B", "B", "B", "B",
          "B", "B", "B", "B", "B", "B", "B", "Mono & DC", "NK", "Platelets",
          "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T", "B", "B",
          "B", "B", "B", "B", "B", "B", "B", "B", "B", "Mono & DC", "NK",
          "Platelets", "T", "T", "T", "T", "T", "T", "T", "T", "T", "T",
          "T", "T"), sample_component = c("S1_128b348aca07cb57", "S1_218e57362e5c1a3d",
          "S1_6d7587035b64afd5", "S1_a0525046264c9183", "S1_a9fec42b5a90b80a",
          "S1_c716219b3865936c", "S1_d375cd40d5330b42", "S1_dc8d14e01732603e",
          "S1_e2055911bf1693bd", "S1_e5bd60fc7282a0f8", "S1_f0dceb55e0820e2d",
          "S1_ff94b67d225a95b9", NA, NA, NA, "S1_2ee5115016f6d3bf", "S1_500dadc305f7fd2d",
          "S1_7e9a3f5ef0362b6d", "S1_844fa3a723a26b8a", "S1_887179ebbb0ff0a2",
          "S1_99aa508551f33cd3", "S1_9f31266d8f62d933", "S1_b720a4c4f4bfac3e",
          "S1_c57ec812b0a4d6e9", "S1_c8d43c718f2181fc", "S1_d0f9201e37a9e091",
          "S2_2ee5115016f6d3bf", "S2_a0525046264c9183", "S2_a9fec42b5a90b80a",
          "S2_c57ec812b0a4d6e9", "S2_c716219b3865936c", "S2_c8d43c718f2181fc",
          "S2_d0f9201e37a9e091", "S2_d375cd40d5330b42", "S2_dc8d14e01732603e",
          "S2_e2055911bf1693bd", "S2_ff94b67d225a95b9", NA, NA, NA, "S2_128b348aca07cb57",
          "S2_218e57362e5c1a3d", "S2_500dadc305f7fd2d", "S2_6d7587035b64afd5",
          "S2_7e9a3f5ef0362b6d", "S2_844fa3a723a26b8a", "S2_887179ebbb0ff0a2",
          "S2_99aa508551f33cd3", "S2_9f31266d8f62d933", "S2_b720a4c4f4bfac3e",
          "S2_e5bd60fc7282a0f8", "S2_f0dceb55e0820e2d"), condition = c("good",
          "good", "good", "good", "good", "good", "good", "good", "good",
          "good", "good", "good", NA, NA, NA, "good", "good", "good", "good",
          "good", "good", "good", "good", "good", "good", "good", "good",
          "good", "good", "good", "good", "good", "good", "good", "good",
          "good", "good", NA, NA, NA, "good", "good", "good", "good", "good",
          "good", "good", "good", "good", "good", "good", "good"), seurat_clusters = c("1",
          "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", NA, NA,
          NA, "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1",
          "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", NA, NA, NA,
          "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1"),
              marker_1 = structure(c(1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L,
              1L, 1L, 1L, NA, NA, NA, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L,
              1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, NA, NA,
              NA, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L), levels = c("CD11b",
              "B2M", "HLA-ABC"), class = "factor"), marker_2 = structure(c(1L,
              1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, NA, NA, NA, 1L,
              1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L,
              1L, 1L, 1L, 1L, 1L, 1L, NA, NA, NA, 1L, 1L, 1L, 1L, 1L, 1L,
              1L, 1L, 1L, 1L, 1L, 1L), levels = c("CD11b", "B2M", "HLA-ABC"
              ), class = "factor"), join_count = c(NA, NA, NA, NA, NA,
              NA, NA, 0, 0, NA, NA, NA, NA, NA, NA, 0, NA, NA, NA, NA,
              0, NA, NA, NA, NA, NA, 0, NA, NA, NA, NA, NA, NA, NA, 0,
              0, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 0, NA, NA,
              NA, NA), join_count_expected_mean = c(NA, NA, NA, NA, NA,
              NA, NA, 0.12, 0.2, NA, NA, NA, NA, NA, NA, 0, NA, NA, NA,
              NA, 0.04, NA, NA, NA, NA, NA, 0, NA, NA, NA, NA, NA, NA,
              NA, 0.12, 0.2, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
              0.04, NA, NA, NA, NA), join_count_expected_sd = c(NA, NA,
              NA, NA, NA, NA, NA, 0.32659863237109, 0.449466574975495,
              NA, NA, NA, NA, NA, NA, 0, NA, NA, NA, NA, 0.196946385566932,
              NA, NA, NA, NA, NA, 0, NA, NA, NA, NA, NA, NA, NA, 0.32659863237109,
              0.449466574975495, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
              NA, 0.196946385566932, NA, NA, NA, NA), join_count_z = c(0,
              0, 0, 0, 0, 0, 0, -0.12, -0.2, 0, 0, 0, NA, NA, NA, 0, 0,
              0, 0, 0, -0.04, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.12,
              -0.2, 0, NA, NA, NA, 0, 0, 0, 0, 0, 0, 0, -0.04, 0, 0, 0,
              0), join_count_p = c(NA, NA, NA, NA, NA, NA, NA, 0.452241573979416,
              0.420740290560897, NA, NA, NA, NA, NA, NA, 0.5, NA, NA, NA,
              NA, 0.484046563147169, NA, NA, NA, NA, NA, 0.5, NA, NA, NA,
              NA, NA, NA, NA, 0.452241573979416, 0.420740290560897, NA,
              NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 0.484046563147169,
              NA, NA, NA, NA), log2_ratio = c(0, 0, 0, 0, 0, 0, 0, 0, 0,
              0, 0, 0, NA, NA, NA, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, NA, NA, NA, 0, 0, 0, 0, 0,
              0, 0, 0, 0, 0, 0, 0), count_1 = c(NA, NA, NA, NA, NA, NA,
              NA, 21L, 49L, NA, NA, NA, NA, NA, NA, 12L, NA, NA, NA, NA,
              19L, NA, NA, NA, NA, NA, 12L, NA, NA, NA, NA, NA, NA, NA,
              21L, 49L, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 19L,
              NA, NA, NA, NA), count_2 = c(NA, NA, NA, NA, NA, NA, NA,
              21L, 49L, NA, NA, NA, NA, NA, NA, 12L, NA, NA, NA, NA, 19L,
              NA, NA, NA, NA, NA, 12L, NA, NA, NA, NA, NA, NA, NA, 21L,
              49L, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 19L, NA,
              NA, NA, NA), p1 = c(NA, NA, NA, NA, NA, NA, NA, 0.010989010989011,
              0.016327890703099, NA, NA, NA, NA, NA, NA, 0.0184331797235023,
              NA, NA, NA, NA, 0.00842945874001775, NA, NA, NA, NA, NA,
              0.0184331797235023, NA, NA, NA, NA, NA, NA, NA, 0.010989010989011,
              0.016327890703099, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
              NA, 0.00842945874001775, NA, NA, NA, NA), p2 = c(NA, NA,
              NA, NA, NA, NA, NA, 0.010989010989011, 0.016327890703099,
              NA, NA, NA, NA, NA, NA, 0.0184331797235023, NA, NA, NA, NA,
              0.00842945874001775, NA, NA, NA, NA, NA, 0.0184331797235023,
              NA, NA, NA, NA, NA, NA, NA, 0.010989010989011, 0.016327890703099,
              NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 0.00842945874001775,
              NA, NA, NA, NA)), row.names = c(NA, -52L), class = c("tbl_df",
          "tbl", "data.frame")),
        hashing =
          structure(list(sample_alias = structure(c(1L, 1L, 1L, 1L, 1L,
          3L, 3L, 3L, 3L, 3L, 3L, 4L, 4L, 4L, 4L, 4L, 2L, 2L, 2L, 2L, 2L
          ), levels = c("S1", "S2", "S11", "S12"), class = "factor"), celltype = c("B",
          "Mono & DC", "NK", "Platelets", "T", "B", "B", "Mono & DC", "NK",
          "Platelets", "T", "B", "Mono & DC", "NK", "Platelets", "T", "B",
          "Mono & DC", "NK", "Platelets", "T"), sample_component = c("S1_d2146defe08567d3",
          NA, NA, NA, "S1_68189b2c75de4098", "S11_19b04397ed7f04ba", "S11_fe556695f452a4bb",
          NA, NA, NA, NA, "S12_63fab986007319b4", NA, NA, NA, "S12_3d23c6539cbead8d",
          "S2_5bfdf506169806f2", NA, NA, NA, "S2_85f5bcdc05a7b286"), condition = c("good",
          NA, NA, NA, "good", "good", "good", NA, NA, NA, NA, "good", NA,
          NA, NA, "good", "good", NA, NA, NA, "good"), seurat_clusters = c("1",
          NA, NA, NA, "1", "1", "1", NA, NA, NA, NA, "1", NA, NA, NA, "1",
          "1", NA, NA, NA, "1"), marker_1 = structure(c(1L, NA, NA, NA,
          1L, 1L, 1L, NA, NA, NA, NA, 1L, NA, NA, NA, 1L, 1L, NA, NA, NA,
          1L), levels = c("CD45", "B2M", "HLA-DR-DP-DQ"), class = "factor"),
              marker_2 = structure(c(1L, NA, NA, NA, 1L, 1L, 1L, NA, NA,
              NA, NA, 1L, NA, NA, NA, 1L, 1L, NA, NA, NA, 1L), levels = c("CD45",
              "B2M", "HLA-DR-DP-DQ"), class = "factor"), join_count = c(153,
              NA, NA, NA, 345, 258, 214, NA, NA, NA, NA, 17, NA, NA, NA,
              NA, 1, NA, NA, NA, 18), join_count_expected_mean = c(143.82,
              NA, NA, NA, 325.28, 245.94, 162.28, NA, NA, NA, NA, 15.12,
              NA, NA, NA, NA, 1.08, NA, NA, NA, 19.49), join_count_expected_sd = c(12.001161559944,
              NA, NA, NA, 19.0899412934105, 16.6192173052278, 13.1978571875678,
              NA, NA, NA, NA, 3.47365263465745, NA, NA, NA, NA, 0.849004170076968,
              NA, NA, NA, 4.54938113229357), join_count_z = c(0.764925957720613,
              NA, NA, NA, 1.03300474825489, 0.725665943137187, 3.91881797665758,
              NA, NA, NA, NA, 0.541217040887393, NA, NA, NA, 0, -0.0800000000000001,
              NA, NA, NA, -0.327517074668311), join_count_p = c(0.222157817923894,
              NA, NA, NA, 0.150800838294423, 0.234021792365565, 4.44921405915763e-05,
              NA, NA, NA, NA, 0.294178996582506, NA, NA, NA, NA, 0.468118627986013,
              NA, NA, NA, 0.371638415373185), log2_ratio = c(0.0892673380970873,
              NA, NA, NA, 0.0849142415955236, 0.0690646698420441, 0.399125588969599,
              NA, NA, NA, NA, 0.169076606803992, NA, NA, NA, 0, -0.111031312388744,
              NA, NA, NA, -0.114737184040853), count_1 = c(934L, NA, NA,
              NA, 1274L, 1582L, 906L, NA, NA, NA, NA, 345L, NA, NA, NA,
              NA, 77L, NA, NA, NA, 409L), count_2 = c(934L, NA, NA, NA,
              1274L, 1582L, 906L, NA, NA, NA, NA, 345L, NA, NA, NA, NA,
              77L, NA, NA, NA, 409L), p1 = c(0.365271802894016, NA, NA,
              NA, 0.594493700419972, 0.366373320981936, 0.46749226006192,
              NA, NA, NA, NA, 0.109211775878443, NA, NA, NA, NA, 0.0275985663082437,
              NA, NA, NA, 0.116690442225392), p2 = c(0.365271802894016,
              NA, NA, NA, 0.594493700419972, 0.366373320981936, 0.46749226006192,
              NA, NA, NA, NA, 0.109211775878443, NA, NA, NA, NA, 0.0275985663082437,
              NA, NA, NA, 0.116690442225392)), row.names = c(NA, -21L), class = c("tbl_df",
          "tbl", "data.frame"))
      )
    )

    # component_proximity_selected
    set.seed(37)
    selected_es_data <- proximity_es_data
    selected_es_data$proximity <- proximity_scores %>%
      filter(as.character(marker_1) == as.character(marker_2))
    expect_no_error(
      component <- component_proximity_selected(
        selected_es_data,
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
      component <- component_dimred_plots(es_data,
        sample_palette = c("red", "black", "blue")
      )
    )

    expect_s3_class(component$PCA, "ggplot")

    # component_proximity_heatmap_sample
    proximity_es_data$pxl_data_processed <- temp

    expect_no_error(
      component <- component_proximity_heatmap_sample(
        proximity_es_data,
        heatmap_gradient = c("red", "cyan"),
        n_markers = 2,
        plot_markers = NULL,
        test_mode = TRUE
      )
    )
    # component_proximity_heatmap_celltype
    expect_no_error(
      component <- component_proximity_heatmap_celltype(
        proximity_es_data,
        heatmap_gradient = c("red", "cyan"),
        n_markers = 2,
        plot_markers = NULL,
        test_mode = TRUE
      )
    )

    # component_hashing
    if (data_type == "hashing") {
      expect_no_error(
        component <- component_hashing(es_data, colors = cherry_gradient)
      )
      expect_s3_class(component$plot, "ggplot")
      expect_no_error(
        ggplot2::ggplot_build(component$plot)
      )

      for (plot in component$sample_confidence_plots) {
        expect_s3_class(plot, "ggplot")
        expect_no_error(
          ggplot2::ggplot_build(plot)
        )
      }
      for (plot in component$heatmap_plots_hash_purity) {
        expect_s3_class(plot, "ggplot")
        expect_no_error(
          ggplot2::ggplot_build(plot)
        )
      }
      for (plot in component$heatmap_plots_hash_fraction) {
        expect_s3_class(plot, "ggplot")
        expect_no_error(
          ggplot2::ggplot_build(plot)
        )
      }
      expect_s3_class(component$table, "datatables")

      qc_metrics_tables_missing_undetermined <-
        qc_metrics_tables

      qc_metrics_tables_missing_undetermined$sample_hash_stats$component_sample_confidence <-
        qc_metrics_tables_missing_undetermined$sample_hash_stats$component_sample_confidence %>%
        filter(sample_alias != "undetermined")
      hashing_es_data <- es_data
      hashing_es_data$qc <- qc_metrics_tables_missing_undetermined
      expect_no_error(
        component <- component_hashing(hashing_es_data,
          colors = cherry_gradient
        )
      )
    }
  })
}
