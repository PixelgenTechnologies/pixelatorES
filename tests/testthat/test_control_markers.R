.test_control_object <- function(markers, control, assay = "PNA") {
  counts <- matrix(
    1,
    nrow = length(markers),
    ncol = 2,
    dimnames = list(markers, c("cell_1", "cell_2"))
  ) %>%
    as("dgCMatrix")
  pna_assay <- SeuratObject::CreateAssay5Object(counts = counts)
  if (!is.null(control)) {
    pna_assay[["control"]] <- set_names(control, markers)
  }
  object <- SeuratObject::CreateSeuratObject(pna_assay, assay = assay)

  return(object)
}

test_that("Control markers are detected from the panel metadata", {
  object <- .test_control_object(
    markers = c("CD3", "mIgG1", "CD4", "rIgG1"),
    control = c("no", "yes", "no", "yes")
  )
  expect_equal(get_control_markers(object), c("mIgG1", "rIgG1"))

  logical_control <- .test_control_object(
    markers = c("CD3", "mIgG1", "CD4", "rIgG1"),
    control = c(FALSE, TRUE, FALSE, TRUE)
  )
  expect_equal(get_control_markers(logical_control), c("mIgG1", "rIgG1"))

  no_controls <- .test_control_object(
    markers = c("CD3", "CD4"),
    control = c("no", "no")
  )
  expect_equal(get_control_markers(no_controls), character())

  no_column <- .test_control_object(markers = c("CD3", "CD4"), control = NULL)
  expect_error(get_control_markers(no_column), "control")
  expect_error(get_control_markers("not a Seurat object"))

  expect_equal(
    pixelatorES:::.get_control_markers_union(
      list(
        S1 = object,
        S2 = .test_control_object(
          markers = c("CD3", "mIgG2a", "rIgG1"),
          control = c("no", "yes", "yes")
        )
      )
    ),
    c("mIgG1", "mIgG2a", "rIgG1")
  )

  expect_equal(
    get_control_markers(suppressMessages(get_test_data(type = "default"))),
    c("mIgG1", "mIgG2a", "mIgG2b")
  )
})

test_that("Downsampling keeps detected control markers", {
  object <- .test_control_object(
    markers = c("CD3", "CD4", "CD8", "mIgG1", "rIgG1"),
    control = c("no", "no", "no", "yes", "yes")
  )
  object$sample_alias <- "S1"

  downsampled <- downsample_data(object, n_cells = 2, n_markers = 3)
  expect_equal(nrow(downsampled), 3)
  expect_equal(
    intersect(rownames(downsampled), c("mIgG1", "rIgG1")),
    c("mIgG1", "rIgG1")
  )
})
