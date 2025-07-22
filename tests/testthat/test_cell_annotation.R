test_that("Reading reference data works as expected", {
  expect_no_error(ref_data <- read_annotation_reference())

  expect_equal(
    LayerData(ref_data, "counts")[1:10, 1:10],
    new("dgCMatrix",
      i = c(
        0L, 1L, 2L, 3L, 4L, 5L, 6L, 7L, 8L, 9L,
        0L, 1L, 2L, 3L, 4L, 5L, 6L, 8L, 9L, 0L, 1L, 2L, 3L, 4L, 5L, 6L,
        7L, 8L, 9L, 0L, 1L, 2L, 3L, 4L, 5L, 6L, 7L, 8L, 9L, 0L, 1L, 2L,
        4L, 5L, 6L, 7L, 8L, 9L, 0L, 1L, 2L, 3L, 4L, 5L, 6L, 7L, 8L, 9L,
        0L, 1L, 2L, 3L, 4L, 5L, 6L, 7L, 9L, 0L, 1L, 2L, 3L, 4L, 5L, 6L,
        7L, 8L, 9L, 0L, 1L, 2L, 3L, 4L, 5L, 6L, 7L, 8L, 9L, 0L, 1L, 2L,
        3L, 4L, 5L, 6L, 7L, 8L, 9L
      ),
      p = c(
        0L, 10L, 19L, 29L, 39L, 48L,
        58L, 67L, 77L, 87L, 97L
      ),
      Dim = c(10L, 10L),
      Dimnames = list(
        c(
          "HLA-ABC", "B2M", "CD11b", "CD11c", "CD18", "CD82", "CD8",
          "TCRab", "HLA-DR", "CD45"
        ),
        c(
          "Resting_PBMC_82c6f1502ceb578c",
          "Resting_PBMC_9524c5af1f13409e", "Resting_PBMC_1a72fb4c1ebd3edf",
          "Resting_PBMC_791355c983aed829", "Resting_PBMC_81415c5542195eff",
          "Resting_PBMC_5ba613cb25de4e5c", "Resting_PBMC_5e08b1ef86b8a4c9",
          "Resting_PBMC_ef963c459c4a8cd6", "Resting_PBMC_89d045aeb08b3401",
          "Resting_PBMC_ce389cbfa22fcb90"
        )
      ),
      x = c(
        1294, 2002, 5, 3,
        132, 1, 11, 7, 2, 1553, 2191, 4690, 3, 6, 295, 5, 4, 6, 1931,
        1368, 3990, 1, 2, 310, 1717, 6, 72, 5, 2098, 17911, 21774, 3,
        12, 18, 1661, 29, 7, 12, 124, 1368, 2357, 2, 146, 25, 1, 10,
        2, 2464, 1482, 2341, 10, 15, 180, 35, 23, 193, 11, 2306, 1252,
        2858, 1, 2, 108, 27, 3390, 97, 1377, 942, 1767, 6, 15, 95, 15,
        477, 6, 1, 976, 2046, 4187, 17, 20, 376, 271, 2872, 95, 25, 2447,
        9771, 11606, 5, 5, 3, 138, 2, 2, 7, 66
      ),
      factors = list()
    )
  )

  expect_equal(
    ref_data$celltype.l1[1:10],
    c(
      Resting_PBMC_82c6f1502ceb578c = "NK", Resting_PBMC_9524c5af1f13409e = "NK",
      Resting_PBMC_1a72fb4c1ebd3edf = "CD4 T", Resting_PBMC_791355c983aed829 = "Platelets",
      Resting_PBMC_81415c5542195eff = "CD4 T", Resting_PBMC_5ba613cb25de4e5c = "CD4 T",
      Resting_PBMC_5e08b1ef86b8a4c9 = "CD8 T", Resting_PBMC_ef963c459c4a8cd6 = "NK",
      Resting_PBMC_89d045aeb08b3401 = "CD8 T", Resting_PBMC_ce389cbfa22fcb90 = "Platelets"
    )
  )
})
