pg_data <- get_test_data()

test_that("ANOVAs work as expected", {
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
})
