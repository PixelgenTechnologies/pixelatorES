pg_data <- get_test_data()

test_that("process_data works as expected", {
  expect_no_error(
    dat_processed <-
      process_data(
        pg_data,
        params = list(
          test_mode = FALSE,
          norm_method = "CLR",
          do_harmonize = FALSE,
          harmonization_vars = "sample_alias",
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
          do_harmonize = FALSE,
          harmonization_vars = "sample_alias",
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
    structure(list(term = c(
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
    )), row.names = c(NA, -4L), class = c(
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
    structure(list(marker = c(
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
    )), row.names = c(
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
    )
  con$close()

  prox_data <-
    prox_data %>%
    mutate(sample_component = paste(component, "1", sep = "_")) %>%
    select(-component) %>%
    left_join(
      FetchData(pg_data, c("condition", "seurat_clusters")) %>%
        as_tibble(rownames = "sample_component"),
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
    structure(
      list(
        marker_1 = c(
          "CD2", "CD2", "CD2", "CD2", "CD2",
          "CD2", "CD2", "CD2", "CD2"
        ),
        marker_2 = c(
          "CD2", "CD2", "CD2",
          "CD4", "CD4", "CD4", "CD3e", "CD3e", "CD3e"
        ),
        term = c(
          "seurat_clusters",
          "condition", "Residuals", "seurat_clusters", "condition", "Residuals",
          "seurat_clusters", "condition", "Residuals"
        ),
        df = c(
          2, 2, NA,
          2, 2, NA, 2, 2, NA
        ),
        sumsq = c(
          0.108108743652551, 0.246958766000748,
          0, 0.0639945478295301, 0.0193967131662747, 0, 0.0261977319993502,
          0.151596778712014, 0
        ),
        meansq = c(
          0.0540543718262753, 0.123479383000374,
          0, 0.031997273914765, 0.00969835658313735, 0, 0.0130988659996751,
          0.0757983893560068, 0
        ),
        statistic = c(
          NA, NA, NA, NA, NA, NA,
          NA, NA, NA
        ),
        p = c(NA, NA, NA, NA, NA, NA, NA, NA, NA),
        eta_squared = c(
          0.304473771081201,
          0.695526228918799, 0, 0.767401128911452, 0.232598871088548, 0,
          0.147348373661998, 0.852651626338002, 0
        ),
        omega_squared = c(
          0.304473771081201,
          0.695526228918799, NA, 0.767401128911452, 0.232598871088548,
          NA, 0.147348373661998, 0.852651626338002, NA
        ),
        p_adj = c(
          NA_real_,
          NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_,
          NA_real_
        )
      ),
      row.names = c(NA, -9L), class = c(
        "tbl_df", "tbl",
        "data.frame"
      )
    )
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
    structure(
      list(
        marker_1 = c(
          "CD2", "CD2", "CD2", "CD2", "CD2",
          "CD2"
        ),
        marker_2 = c("CD2", "CD2", "CD4", "CD4", "CD3e", "CD3e"),
        term = c(
          "seurat_clusters", "Residuals", "seurat_clusters",
          "Residuals", "seurat_clusters", "Residuals"
        ),
        df = c(
          2, 2, 2,
          2, 2, 2
        ),
        sumsq = c(
          0.108108743652551, 0.246958766000748, 0.0639945478295301,
          0.0193967131662747, 0.0261977319993502, 0.151596778712014
        ),
        meansq = c(
          0.0540543718262753,
          0.123479383000374, 0.031997273914765, 0.00969835658313735, 0.0130988659996751,
          0.0757983893560068
        ),
        statistic = c(
          0.437760300649637, NA, 3.29924700545648,
          NA, 0.172811930582758, NA
        ),
        p = c(
          0.695526228918799, NA, 0.232598871088548,
          NA, 0.852651626338002, NA
        ),
        eta_squared = c(
          0.304473771081201,
          0.695526228918799, 0.767401128911452, 0.232598871088548, 0.147348373661998,
          0.852651626338002
        ),
        omega_squared = c(
          -0.290149250741628, NA,
          0.479084948710154, NA, -0.494489580265652, NA
        ),
        p_adj = c(
          1,
          NA, 0.697796613265645, NA, 1, NA
        )
      ),
      row.names = c(NA, -6L),
      class = c(
        "tbl_df",
        "tbl", "data.frame"
      )
    )
  )
})
