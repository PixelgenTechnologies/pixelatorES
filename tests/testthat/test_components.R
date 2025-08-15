library(dplyr)
pg_data <- get_test_data()
pg_data_small <- get_test_data(concatenate = FALSE)
sample_qc_metrics <- get_test_qc_metrics()
sample_sheet <- read_samplesheet(test_samplesheet())[1, ]
qc_metrics_tables <-
  get_qc_metrics(pg_data, sample_qc_metrics[1], sample_sheet)

test_that("Components work as expected", {
  # component_control_markers
  expect_no_error(
    component <- component_control_markers(pg_data)
  )

  expect_s3_class(component$p1, "ggplot")
  expect_s3_class(component$p2, "ggplot")
  expect_s3_class(component$tabl, "datatables")

  # component_cell_recovery
  expect_no_error(component <- component_cell_recovery(sample_qc_metrics, sample_levels = NULL))

  # component_node_degree
  expect_no_error(
    component <- component_node_degree(pg_data, sample_levels = NULL)
  )

  expect_s3_class(component$plot, "ggplot")
  expect_s3_class(component$table, "datatables")

  # component_node_edge_count
  expect_no_error(
    component <- component_node_edge_count(sample_qc_metrics, sample_levels = NULL)
  )

  for (plot in component$plots) expect_s3_class(plot, "ggplot")
  expect_s3_class(component$table, "datatables")

  # component_sequencing_saturation
  expect_no_error(
    component <- component_sequencing_saturation(qc_metrics_tables, sample_levels = NULL)
  )

  for (plot in component$plots) expect_s3_class(plot, "ggplot")
  expect_s3_class(component$table, "datatables")

  # component_sequencing_reads_per_cell
  expect_no_error(component <- component_sequencing_reads_per_cell(pg_data))
  expect_s3_class(component$plot, "ggplot")
  expect_s3_class(component$table, "datatables")

  # component_abundance_per_celltype
  temp <-
    pg_data %>%
    subset(features = rownames(pg_data)[1:3])

  temp[["l1_annotation_summary"]] <- sample(
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
      sample_palette = c("red", "black")
    )
  )

  expect_s3_class(component[[1]], "ggplot")

  # component_proximity_per_marker
  temp <-
    pg_data_small %>%
    subset(features = rownames(pg_data_small)[1:3])

  temp[["l1_annotation_summary"]] <- sample(
    c("CD4 T", "CD8 T", "B"),
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
      proximity_score = "log2_ratio",
      sample_palette = c("red", "black")
    )
  )
  expect_s3_class(component[[1]], "ggplot")
  expect_equal(
    component[[1]]$data,
    structure(
      list(
        sample_alias = c("S1", "S1", "S1"),
        marker_1 = structure(c(1L,
                               1L, 1L),
                             levels = c("CD11b", "B2M", "HLA-ABC"), class = "factor"),
        marker_2 = structure(c(1L, 1L, 1L),
                             levels = c("CD11b", "B2M",
                                        "HLA-ABC"),
                             class = "factor"),
        join_count = c(0, 57, 0),
        join_count_expected_mean = c(0.03, 45.18, 0),
        join_count_expected_sd = c(0.171446607997765,
                                   6.58338895258438, 0),
        join_count_z = c(-0.03, 1.79542786931341,
                         0),
        join_count_p = c(0.488033526585887, 0.0362927777357692,
                         0.5),
        log2_ratio = c(0, 0.335277648546382, 0),
        sample_component = c("c3c393e9a17c1981",
                             "0a45497c6bfbfb22", "2708240b908e2eba"),
        count_1 = c(17L,
                    929L, 12L),
        count_2 = c(17L, 929L, 12L),
        p1 = c(0.0021783700666325,
               0.312163978494624, 0.00216723857684667),
        p2 = c(0.0021783700666325,
               0.312163978494624, 0.00216723857684667),
        condition = c("good",
                      "good", "good"),
        seurat_clusters = c("1", "1", "1"),
        l1_annotation_summary = c("CD4 T",
                                  "CD4 T", "B")),
      class = c("tbl_df", "tbl", "data.frame"),
      row.names = c(NA,
                    -3L))
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

  expect_named(component, expected = c("B2M/B2M", "CD11b/CD11b", "HLA-ABC/HLA-ABC"))
  expect_s3_class(component[[1]], "ggplot")
  expect_equal(component[[1]]$data,
               structure(list(sample_alias = c("S1", "S1", "S1", "S1", "S1",
                                               "S1", "S1"),
                              l1_annotation_summary = c("B", "CD4 T", "CD4 T",
                                                        "CD4 T", "CD8 T", "Mono", "NK"),
                              contrast = c("B2M/B2M", "B2M/B2M",
                                           "B2M/B2M", "B2M/B2M", "B2M/B2M", NA, NA),
                              marker_1 = structure(c(2L,
                                                     2L, 2L, 2L, 2L, NA, NA), levels = c("CD11b", "B2M", "HLA-ABC"
                                                     ), class = "factor"),
                              marker_2 = structure(c(2L, 2L, 2L, 2L,
                                                     2L, NA, NA),
                                                   levels = c("CD11b", "B2M", "HLA-ABC"), class = "factor"),
                              join_count = c(710, 1572, 1653, 216, 3311, NA, NA),
                              join_count_expected_mean = c(675.83,
                                                           1296.95, 1614.79, 73.17, 3175.65, NA, NA),
                              join_count_expected_sd = c(28.2446222645779,
                                                         38.9973775342068, 40.0379352942464,
                                                         8.01772657266006, 57.9514556514435,
                                                         NA, NA),
                              join_count_z = c(1.20978782013499, 7.05303836799637,
                                               0.954344916119861, 17.8142767411203, 2.33557549984732, NA,
                                               NA),
                              join_count_p = c(0.113180160531595, 8.75262779783333e-13,
                                               0.169954539545749, 2.73805633201208e-71, 0.00975668836471254,
                                               NA, NA),
                              log2_ratio = c(0.0711586316613331, 0.277478355695801,
                                             0.0337401669294608, 1.56170714842085, 0.0602150941888999,
                                             NA, NA),
                              sample_component = c("2708240b908e2eba", "c3c393e9a17c1981",
                                                   "efe0ed189cb499fc", "0a45497c6bfbfb22", "d4074c845bb62800",
                                                   NA, NA),
                              count_1 = c(3448L, 5307L, 6082L, 1182L, 9753L, NA,
                                          NA),
                              count_2 = c(3448L, 5307L, 6082L, 1182L, 9753L, NA, NA
                              ),
                              p1 = c(0.622719884413943, 0.680035879036392, 0.733301181577044,
                                     0.397177419354839, 0.764761232651141, NA, NA),
                              p2 = c(0.622719884413943,
                                     0.680035879036392, 0.733301181577044, 0.397177419354839,
                                     0.764761232651141, NA, NA),
                              condition = c("good", "good",
                                            "good", "good", "good", NA, NA),
                              seurat_clusters = c("1",
                                                  "1", "1", "1", "1", NA, NA)),
                         class = c("tbl_df", "tbl",
                                   "data.frame"),
                         row.names = c(NA, -7L)))

  # component_dimred_plots
  expect_no_error(
    component <- component_dimred_plots(pg_data, sample_palette = c("red", "black", "blue"))
  )

  expect_s3_class(component$PCA, "ggplot")
})
