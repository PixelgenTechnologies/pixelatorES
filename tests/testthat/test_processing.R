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
    structure(
      list(
        term = c(
          "seurat_clusters", "condition",
          "seurat_clusters:condition",
          "Residuals"
        ),
        df = c(2, 2, 4, 4731),
        sumsq = c(
          4.99208396221759, 31.1997967847596,
          250.687103881314, 14465.9717548549
        ),
        meansq = c(
          2.4960419811088, 15.5998983923798,
          62.6717759703285, 3.05769853199215
        ),
        statistic = c(
          0.816313954758181, 5.10184317687334,
          20.4963881542294, NA
        ),
        p = c(
          0.442120355858416, 0.00611904451531066,
          9.23535811569634e-17, NA
        ),
        eta_squared = c(
          0.000338380971269318, 0.00211483172545489,
          0.0169924517171721, 0.980554335586104
        ),
        omega_squared = c(
          -7.61263263787104e-05, 0.00169995631418742,
          0.0161600562076557, NA
        )
      ),
      row.names = c(NA, -4L),
      class = c("tbl_df", "tbl", "data.frame")
    )
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
    structure(
      list(
        term = c("seurat_clusters", "condition", "Residuals"),
        df = c(2, 2, 4735),
        sumsq = c(
          4.99208396221759, 31.1997967847596,
          14716.6588587362
        ),
        meansq = c(
          2.4960419811088, 15.5998983923798,
          3.1080588930805
        ),
        statistic = c(
          0.80308709293307, 5.01917721929593,
          NA
        ),
        p = c(0.44800498702214, 0.00664517456335452, NA),
        eta_squared = c(
          0.000338380971269319,
          0.00211483172545489, 0.997546787303276
        ),
        omega_squared = c(
          -8.29518325896995e-05,
          0.00169312474641416, NA
        )
      ),
      row.names = c(NA, -3L),
      class = c(
        "tbl_df",
        "tbl", "data.frame"
      )
    )
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
    structure(
      list(
        marker = c(
          "B2M", "B2M", "B2M", "B2M", "CD11b",
          "CD11b", "CD11b", "CD11b", "HLA-ABC",
          "HLA-ABC", "HLA-ABC", "HLA-ABC"
        ), term = c(
          "seurat_clusters", "condition", "seurat_clusters:condition",
          "Residuals", "seurat_clusters", "condition", "seurat_clusters:condition",
          "Residuals", "seurat_clusters", "condition", "seurat_clusters:condition",
          "Residuals"
        ),
        df = c(2, 2, 4, 21, 2, 2, 4, 21, 2, 2, 4, 21),
        sumsq = c(
          0.067520955412666, 1.12116220280441, 1.43807030587216,
          8.25850111822215, 0.742152037086004, 13.3852376955609, 14.956850510018,
          76.055635767835, 0.0170435273467141, 0.524459985602968, 0.478870019773764,
          3.01002385930892
        ),
        meansq = c(
          0.033760477706333, 0.560581101402204,
          0.35951757646804, 0.393261958010579, 0.371076018543002, 6.69261884778044,
          3.73921262750449, 3.62169694132547, 0.00852176367335704,
          0.262229992801484, 0.119717504943441, 0.143334469490901
        ),
        statistic = c(
          0.085847301063951, 1.42546485868619, 0.914193628813619,
          NA, 0.102459157835331, 1.8479234889629, 1.03244768628708,
          NA, 0.059453693892508, 1.8294970758456, 0.835231785966466,
          NA
        ),
        p = c(
          0.918054739230943, 0.262722223998828, 0.473958893685145,
          NA, 0.90306342744124, 0.182285186612611, 0.413784459166228,
          NA, 0.942437188335321, 0.1851660404491, 0.518039937844055,
          NA
        ),
        p_adj = c(1, 1, 1, NA, 1, 1, 1, NA, 1, 1, 1, NA),
        eta_squared = c(
          0.00620297439091485,
          0.10299825275803, 0.132111775153981, 0.758686997697075, 0.00705871135906503,
          0.1273088594305, 0.142256687734008, 0.723375741476428, 0.00422874612325008,
          0.13012612270933, 0.118814591514086, 0.746830539653333
        ),
        omega_squared = c(
          -0.0637497811026809, 0.029670416812963,
          -0.0119676666419378, NA, -0.0597751730608438, 0.0564707152187877,
          0.00432195610966636, NA, -0.0646005590633903, 0.0569732446909948,
          -0.0226339068546109, NA
        )
      ),
      row.names = c(NA, -12L),
      class = c("tbl_df", "tbl", "data.frame")
    )
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
    structure(
      list(
        marker = c(
          "B2M", "B2M", "CD11b", "CD11b", "HLA-ABC",
          "HLA-ABC"
        ),
        term = c(
          "seurat_clusters", "Residuals", "seurat_clusters",
          "Residuals", "seurat_clusters", "Residuals"
        ),
        df = c(
          2, 27, 2,
          27, 2, 27
        ),
        sumsq = c(
          0.067520955412666, 10.8177336268987, 0.742152037086004,
          104.397723973414, 0.0170435273467141, 4.01335386468565
        ),
        meansq = c(
          0.033760477706333, 0.400656800996249, 0.371076018543002,
          3.8665823693857, 0.00852176367335704, 0.148642735729098
        ),
        statistic = c(
          0.0842628344817466, NA, 0.0959700280746785,
          NA, 0.0573305088308384, NA
        ),
        p = c(
          0.919430380488634, NA, 0.908799741070685,
          NA, 0.944396556962291, NA
        ),
        p_adj = c(1, NA, 1, NA, 1, NA),
        eta_squared = c(
          0.00620297439091485, 0.993797025609085, 0.00705871135906503,
          0.992941288640935, 0.00422874612325008, 0.99577125387675
        ),
        omega_squared = c(
          -0.0650184660908417, NA,
          -0.0641339311962768, NA, -0.0670589263428768, NA
        )
      ),
      row.names = c(
        NA,
        -6L
      ),
      class = c("tbl_df", "tbl", "data.frame")
    )
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
