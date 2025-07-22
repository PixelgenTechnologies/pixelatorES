library(Seurat)

test_that("get_top_markers works as expected", {
  set.seed(37)
  seur <-
    CreateSeuratObject(
      counts = matrix(
        c(
          rpois(1000, 40),
          rpois(1000, 5)
        )[sample(1:2000, 2000)],
        nrow = 100,
        ncol = 20,
        dimnames = list(
          paste0("Feature", 1:100),
          paste0("Cell", 1:20)
        )
      ) %>%
        as("dgCMatrix")
    ) %>%
    AddMetaData(
      metadata = data.frame(
        sample_alias = rep(c("A", "B"), each = 10),
        sample_type = rep(
          c("Unstimulated", "Stimulated"),
          each = 5,
          times = 2
        ),
        row.names = paste0("Cell", 1:20)
      )
    )

  expect_no_error(res <- get_top_markers(seur, group = "sample_alias"))
  expect_equal(
    res,
    structure(
      list(
        sample_alias = c("A", "B"),
        top3_fraction = c(
          0.0447203890239667,
          0.0458593463092178
        ),
        top5_fraction = c(
          0.0730288294546718,
          0.0742287917737789
        ),
        top_markers = c(
          "Feature35, Feature95, Feature69, Feature91, Feature31",
          "Feature94, Feature53, Feature80, Feature92, Feature38"
        )
      ),
      row.names = c(NA, -2L),
      class = c("tbl_df", "tbl", "data.frame")
    )
  )
})
