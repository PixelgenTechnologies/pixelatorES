pg_data <- get_test_data()

test_that("process_data works as expected", {
  expect_no_error(
    dat_processed <-
      process_data(
        pg_data,
        params = list(
          test_mode = FALSE,
          norm_method = "CLR",
          clustering_resolution = 1,
          annotation_method = "nmf"
        )
      )
  )

  # With few cells
  expect_no_error(
    dat_processed <-
      process_data(
        pg_data |>
          subset(cells = 1:5),
        params = list(
          test_mode = FALSE,
          norm_method = "CLR",
          clustering_resolution = 1,
          annotation_method = "nmf"
        )
      )
  )
})

test_that("Abundance ANOVAs work as expected", {
  set.seed(37)
  test_data <-
    LayerData(pg_data) %>%
    as_tibble(rownames = "marker") %>%
    pivot_longer(
      cols = -marker,
      names_to = "comp_id",
      values_to = "x"
    ) %>%
    left_join(
      pg_data[[]] %>%
        as_tibble(rownames = "comp_id"),
      by = "comp_id"
    )

  expect_no_error(
    anova_res <- .run_anova(test_data, c("seurat_clusters", "condition"), response = "x", tidy = TRUE)
  )

  expect_equal(
    anova_res,
    structure(list(
      term = c(
        "seurat_clusters", "condition", "seurat_clusters:condition",
        "Residuals"
      ), df = c(2, 2, 4, 7575), sumsq = c(
        35.9647703956345,
        1.8367089075964, 90.6173653391415, 16325.5725550414
      ), meansq = c(
        17.9823851978173,
        0.918354453798199, 22.6543413347854, 2.1551910963751
      ), statistic = c(
        8.34375440213283,
        0.426112772710881, 10.511523257909, NA
      ), p = c(
        0.000240070638436281,
        0.653058344963662, 1.72042453371239e-08, NA
      ), eta_squared = c(
        0.00218577787735599,
        0.000111626951964476, 0.00550731814171746, 0.992195277028962
      ),
      omega_squared = c(
        0.0019235601742038, -0.000150319108517162,
        0.00498273398947364, NA
      )
    ), row.names = c(NA, -4L), class = c(
      "tbl_df",
      "tbl", "data.frame"
    ))
  )

  # Using formula
  expect_no_error(
    anova_res_form <- .run_anova(test_data, formula_str = "x ~ seurat_clusters*condition", tidy = TRUE)
  )
  expect_equal(anova_res, anova_res_form)

  expect_no_error(
    anova_res_form <- .run_anova(test_data, formula_str = "x ~ seurat_clusters + condition", tidy = TRUE)
  )

  expect_equal(
    anova_res_form,
    structure(list(term = c(
      "seurat_clusters", "condition", "Residuals"
    ), df = c(2, 2, 7579), sumsq = c(
      35.9647703956345, 1.8367089075964,
      16416.1899203806
    ), meansq = c(
      17.9823851978173, 0.918354453798199,
      2.16601001720287
    ), statistic = c(
      8.30207850148323, 0.423984398273531,
      NA
    ), p = c(0.000250263143574241, 0.654449613374602, NA), eta_squared = c(
      0.00218577787735599,
      0.000111626951964476, 0.997702595170679
    ), omega_squared = c(
      0.00192224403143963,
      -0.000151633887832061, NA
    )), row.names = c(NA, -3L), class = c(
      "tbl_df",
      "tbl", "data.frame"
    ))
  )

  expect_no_error(
    anova_res <-
      run_abundance_anova(
        pg_data %>%
          subset(features = rownames(pg_data)[1:3]),
        vars = c("seurat_clusters", "condition"),
        layer = "data",
        mc_cores = 1,
        p_adj_method = "bonferroni"
      )
  )

  expect_equal(
    anova_res,
    structure(list(
      marker = c(
        "B2M", "B2M", "B2M", "B2M", "CD11b",
        "CD11b", "CD11b", "CD11b", "HLA-ABC", "HLA-ABC", "HLA-ABC", "HLA-ABC"
      ), term = c(
        "seurat_clusters", "condition", "seurat_clusters:condition",
        "Residuals", "seurat_clusters", "condition", "seurat_clusters:condition",
        "Residuals", "seurat_clusters", "condition", "seurat_clusters:condition",
        "Residuals"
      ), df = c(2, 2, 4, 39, 2, 2, 4, 39, 2, 2, 4, 39),
      sumsq = c(
        1.91763037970324, 11.7915652723369, 5.93708392939578,
        90.4408208295632, 0.347454577494998, 1.44972321689678, 1.5356078656992,
        23.614001145579, 1.5982026884192, 9.17401760870464, 3.94449591837103,
        77.6506629031678
      ), meansq = c(
        0.958815189851618, 5.89578263616847,
        1.48427098234895, 2.31899540588623, 0.173727288747499, 0.724861608448389,
        0.383901966424799, 0.605487208861001, 0.799101344209602,
        4.58700880435232, 0.986123979592757, 1.99104263854276
      ), statistic = c(
        0.413461444303808,
        2.54238650978065, 0.640049125833309, NA, 0.286921484393208,
        1.19715428805168, 0.634038111468891, NA, 0.401348182475118,
        2.30382248755332, 0.495280191645969, NA
      ), p = c(
        0.664221837192148,
        0.091649892551649, 0.63711007725158, NA, 0.752141284421668,
        0.312908598312421, 0.641239763487647, NA, 0.672150013688813,
        0.113315283014414, 0.739222451697284, NA
      ), p_adj = c(
        1, 0.824849032964841,
        1, NA, 1, 1, 1, NA, 1, 1, 1, NA
      ), eta_squared = c(
        0.0174192105391454,
        0.107111234906854, 0.0539307866882702, 0.82153876786573,
        0.0128941005100426, 0.0537994836769086, 0.0569866781065,
        0.876319737706549, 0.0173026744254162, 0.0993209691152863,
        0.0427044261297444, 0.840671930329553
      ), omega_squared = c(
        -0.024201182438547,
        0.0636404494665305, -0.0297038845614599, NA, -0.0313411459167249,
        0.00866530287296291, -0.0321694307075108, NA, -0.025264121042638,
        0.0550235181442365, -0.0426000621984023, NA
      )
    ), row.names = c(
      NA,
      -12L
    ), class = c("tbl_df", "tbl", "data.frame"))
  )

  # When there is only one sample

  test_data$condition <- "one"

  expect_no_error(
    anova_res <- .run_anova(test_data, c("seurat_clusters", "condition"), response = "x", tidy = TRUE)
  )
  expect_null(anova_res)

  pg_data_single_cond <-
    pg_data %>%
    subset(features = rownames(pg_data)[1:3])

  pg_data_single_cond$condition <- "one"

  expect_no_error(
    anova_res <-
      run_abundance_anova(
        pg_data_single_cond,
        vars = c("seurat_clusters", "condition"),
        layer = "data",
        mc_cores = 1,
        p_adj_method = "bonferroni"
      )
  )

  expect_equal(
    anova_res,
    structure(list(marker = c(
      "B2M", "B2M", "CD11b", "CD11b", "HLA-ABC",
      "HLA-ABC"
    ), term = c(
      "seurat_clusters", "Residuals", "seurat_clusters",
      "Residuals", "seurat_clusters", "Residuals"
    ), df = c(
      2, 45, 2,
      45, 2, 45
    ), sumsq = c(
      1.91763037970324, 108.169470031296, 0.347454577494998,
      26.599332228175, 1.5982026884192, 90.7691764302435
    ), meansq = c(
      0.958815189851618,
      2.40376600069546, 0.173727288747499, 0.591096271737223, 0.799101344209602,
      2.01709280956097
    ), statistic = c(
      0.398880419131567, NA, 0.293906926932422,
      NA, 0.396164886623899, NA
    ), p = c(
      0.673420063340052, NA, 0.746765625646715,
      NA, 0.675219409165547, NA
    ), p_adj = c(1, NA, 1, NA, 1, NA), eta_squared = c(
      0.0174192105391454,
      0.982580789460855, 0.0128941005100426, 0.987105899489958, 0.0173026744254161,
      0.982697325574584
    ), omega_squared = c(
      -0.0256901001287626, NA,
      -0.0303123505765876, NA, -0.0258091493329032, NA
    )), row.names = c(
      NA,
      -6L
    ), class = c("tbl_df", "tbl", "data.frame"))
  )
})

test_that("Proximity ANOVAs work as expected", {
  set.seed(37)

  con <- pixelatorR::PixelDB$new(pixelatorR::FSMap(pg_data)$pxl_file[1])
  prox_data <-
    con$proximity() %>%
    group_by(marker_1, marker_2) %>%
    filter(any(join_count > 0)) %>%
    ungroup() %>%
    filter(
      marker_1 == "CD2",
      marker_2 %in% c("CD3e", "CD2", "CD4")
    ) %>%
    mutate(component = paste0("S1_", component)) %>%
    distinct()
  con$close()

  prox_data <-
    prox_data %>%
    rename(sample_component = component) %>%
    left_join(
      FetchData(pg_data, c("condition", "seurat_clusters")) %>%
        as_tibble(rownames = "sample_component") %>%
        distinct(),
      by = c("sample_component")
    )

  expect_no_error(
    anova_res <-
      run_proximity_anova(
        prox_data,
        vars = c("seurat_clusters", "condition"),
        mc_cores = 1,
        p_adj_method = "bonferroni"
      )
  )

  expect_equal(
    anova_res,
    structure(list(marker_1 = c(
      "CD2", "CD2", "CD2", "CD2", "CD2",
      "CD2", "CD2", "CD2", "CD2", "CD2", "CD2", "CD2"
    ), marker_2 = c(
      "CD2",
      "CD2", "CD2", "CD2", "CD3e", "CD3e", "CD3e", "CD3e", "CD4", "CD4",
      "CD4", "CD4"
    ), term = c(
      "seurat_clusters", "condition", "seurat_clusters:condition",
      "Residuals", "seurat_clusters", "condition", "seurat_clusters:condition",
      "Residuals", "seurat_clusters", "condition", "seurat_clusters:condition",
      "Residuals"
    ), df = c(2, 2, 3, 5, 2, 2, 3, 5, 2, 2, 3, 5), sumsq = c(
      0.187843913016432,
      0.516512171201667, 0.278355527318758, 1.20860760959124, 0.243313643187748,
      0.211955745807452, 0.26994585311777, 2.14445077210472, 0.570444836293754,
      0.238854290366726, 0.166715921445579, 3.00302929595981
    ), meansq = c(
      0.0939219565082161,
      0.258256085600834, 0.0927851757729192, 0.241721521918248, 0.121656821593874,
      0.105977872903726, 0.0899819510392568, 0.428890154420943, 0.285222418146877,
      0.119427145183363, 0.055571973815193, 0.600605859191963
    ), statistic = c(
      0.388554381764902,
      1.06840335751393, 0.383851528968529, NA, 0.283654964656688, 0.247097938274684,
      0.20980185744936, NA, 0.474891168278989, 0.198844455736807, 0.0925265262812418,
      NA
    ), p = c(
      0.696864997396638, 0.410832860473355, 0.769744957795981,
      NA, 0.764383869109485, 0.790068349792672, 0.885579247136356,
      NA, 0.647399434549428, 0.825858006315813, 0.960965237706594,
      NA
    ), eta_squared = c(
      0.0857218387924922, 0.235708319546326, 0.127026461792938,
      0.551543379868244, 0.0847881397982403, 0.0738607715174249, 0.0940687354487701,
      0.747282353235565, 0.143362271683271, 0.0600280544053098, 0.0418984829093975,
      0.754711191002022
    ), omega_squared = c(
      -0.121493703574381, 0.0135916866413699,
      -0.183642234398651, NA, -0.18628352352955, -0.195790076026196,
      -0.308233226346629, NA, -0.13773254595765, -0.210137758410589,
      -0.357036363818167, NA
    ), p_adj = c(
      1, 1, 1, NA, 1, 1, 1, NA,
      1, 1, 1, NA
    )), row.names = c(NA, -12L), class = c(
      "tbl_df", "tbl",
      "data.frame"
    ))
  )

  # Only one sample
  prox_data$condition <- "one"

  expect_no_error(
    anova_res <-
      run_proximity_anova(
        prox_data,
        vars = c("seurat_clusters", "condition"),
        mc_cores = 1,
        p_adj_method = "bonferroni"
      )
  )

  expect_equal(
    anova_res,
    structure(list(marker_1 = c(
      "CD2", "CD2", "CD2", "CD2", "CD2",
      "CD2"
    ), marker_2 = c("CD2", "CD2", "CD3e", "CD3e", "CD4", "CD4"), term = c(
      "seurat_clusters", "Residuals", "seurat_clusters",
      "Residuals", "seurat_clusters", "Residuals"
    ), df = c(
      2, 10, 2,
      10, 2, 10
    ), sumsq = c(
      0.187843913016432, 2.00347530811166, 0.243313643187748,
      2.62635237102994, 0.570444836293754, 3.40859950777212
    ), meansq = c(
      0.0939219565082161,
      0.200347530811166, 0.121656821593874, 0.262635237102994, 0.285222418146877,
      0.340859950777212
    ), statistic = c(
      0.468795178697462, NA, 0.463215914725737,
      NA, 0.836773042701342, NA
    ), p = c(
      0.638839455503186, NA, 0.642108172483754,
      NA, 0.461302692101845, NA
    ), eta_squared = c(
      0.0857218387924921,
      0.914278161207508, 0.0847881397982403, 0.91521186020176, 0.143362271683271,
      0.856637728316729
    ), omega_squared = c(
      -0.0889969927596778, NA,
      -0.0900158728025783, NA, -0.0257586876157197, NA
    ), p_adj = c(
      1,
      NA, 1, NA, 1, NA
    )), row.names = c(NA, -6L), class = c(
      "tbl_df",
      "tbl", "data.frame"
    ))
  )
})

test_that("Proximity scores are zero-padded as expected", {
  proximity_scores <- tibble(
    sample_component = c("c1", "c1", "c2"),
    sample_alias = factor(c("S1", "S1", "S1"), levels = c("S1", "S2")),
    condition = "unstim",
    seurat_clusters = "1",
    celltype = "T",
    marker_1 = factor(c("CD3", "CD3", "CD4"), levels = c("CD3", "CD4")),
    marker_2 = factor(c("CD3", "CD4", "CD4"), levels = c("CD3", "CD4")),
    log2_ratio = c(0.5, 0.2, 0.1),
    join_count_z = c(1.2, 0.3, 0.8)
  )

  completed <- complete_proximity_scores(proximity_scores, only_self = FALSE) %>%
    arrange(sample_component, marker_1, marker_2)

  expect_equal(
    completed,
    structure(list(
      sample_component = c("c1", "c1", "c1", "c2", "c2", "c2"),
      sample_alias = structure(
        c(1L, 1L, 1L, 1L, 1L, 1L),
        levels = c("S1", "S2"),
        class = "factor"
      ),
      condition = c("unstim", "unstim", "unstim", "unstim", "unstim", "unstim"),
      seurat_clusters = c("1", "1", "1", "1", "1", "1"),
      celltype = c("T", "T", "T", "T", "T", "T"),
      marker_1 = structure(
        c(1L, 1L, 2L, 1L, 1L, 2L),
        levels = c("CD3", "CD4"),
        class = "factor"
      ),
      marker_2 = structure(
        c(1L, 2L, 2L, 1L, 2L, 2L),
        levels = c("CD3", "CD4"),
        class = "factor"
      ),
      log2_ratio = c(0.5, 0.2, 0, 0, 0, 0.1),
      join_count_z = c(1.2, 0.3, 0, 0, 0, 0.8)
    ), row.names = c(NA, -6L), class = c("tbl_df", "tbl", "data.frame"))
  )

  completed_self <- complete_proximity_scores(proximity_scores, only_self = TRUE) %>%
    arrange(sample_component, marker_1, marker_2)

  expect_equal(
    completed_self,
    structure(list(
      sample_component = c("c1", "c1", "c2", "c2"),
      sample_alias = structure(
        c(1L, 1L, 1L, 1L),
        levels = c("S1", "S2"),
        class = "factor"
      ),
      condition = c("unstim", "unstim", "unstim", "unstim"),
      seurat_clusters = c("1", "1", "1", "1"),
      celltype = c("T", "T", "T", "T"),
      marker_1 = structure(
        c(1L, 2L, 1L, 2L),
        levels = c("CD3", "CD4"),
        class = "factor"
      ),
      marker_2 = structure(
        c(1L, 2L, 1L, 2L),
        levels = c("CD3", "CD4"),
        class = "factor"
      ),
      log2_ratio = c(0.5, 0, 0, 0.1),
      join_count_z = c(1.2, 0, 0, 0.8)
    ), row.names = c(NA, -4L), class = c("tbl_df", "tbl", "data.frame"))
  )

  component_meta <- tibble(
    sample_component = c("c1", "c2", "c3"),
    sample_alias = factor(c("S1", "S1", "S2"), levels = c("S1", "S2")),
    condition = c("unstim", "unstim", "stim"),
    seurat_clusters = c("1", "1", "2"),
    celltype = c("T", "T", "B")
  )

  completed_all <- complete_proximity_scores(
    proximity_scores,
    only_self = FALSE,
    component_meta = component_meta
  ) %>%
    arrange(sample_component, marker_1, marker_2)

  expect_equal(
    completed_all,
    structure(list(
      sample_component = c(
        "c1", "c1", "c1", "c2", "c2", "c2", "c3", "c3", "c3"
      ),
      sample_alias = structure(
        c(1L, 1L, 1L, 1L, 1L, 1L, 2L, 2L, 2L),
        levels = c("S1", "S2"),
        class = "factor"
      ),
      condition = c(
        "unstim", "unstim", "unstim",
        "unstim", "unstim", "unstim",
        "stim", "stim", "stim"
      ),
      seurat_clusters = c("1", "1", "1", "1", "1", "1", "2", "2", "2"),
      celltype = c("T", "T", "T", "T", "T", "T", "B", "B", "B"),
      marker_1 = structure(
        c(1L, 1L, 2L, 1L, 1L, 2L, 1L, 1L, 2L),
        levels = c("CD3", "CD4"),
        class = "factor"
      ),
      marker_2 = structure(
        c(1L, 2L, 2L, 1L, 2L, 2L, 1L, 2L, 2L),
        levels = c("CD3", "CD4"),
        class = "factor"
      ),
      log2_ratio = c(0.5, 0.2, 0, 0, 0, 0.1, 0, 0, 0),
      join_count_z = c(1.2, 0.3, 0, 0, 0, 0.8, 0, 0, 0)
    ), row.names = c(NA, -9L), class = c("tbl_df", "tbl", "data.frame"))
  )
})
