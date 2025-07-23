test_that("`order_cd_markers` works as expected", {
  markers <- c(
    "CD3", "CD4", "CD8", "CD14", "CD16", "CD19", "CD20", "CD34",
    "Control1", "Control2", "ABC", "WTF"
  )

  expect_equal(
    order_cd_markers(markers, control_markers = c("Control1", "Control2")),
    c(
      "CD3", "CD4", "CD8", "CD14", "CD16", "CD19", "CD20", "CD34", "ABC",
      "WTF", "Control1", "Control2"
    )
  )

  # randomize marker order
  expect_equal(
    order_cd_markers(sample(markers),
      control_markers = c("Control1", "Control2")
    ),
    c(
      "CD3", "CD4", "CD8", "CD14", "CD16", "CD19", "CD20", "CD34", "ABC",
      "WTF", "Control1", "Control2"
    )
  )


  labels <-
    c(
      "PC_1", "pc_1", "PC 1", "PC1", "pc1",
      "UMAP_1", "umap_1", "UMAP 1", "UMAP1", "umap1",
      "hArmony_1", "harmony_1", "hArmony 1", "harmony 1"
    )
  expect_equal(
    beaut_label(labels),
    c(
      "PC 1", "PC 1", "PC 1", "PC1", "PC1", "UMAP 1", "UMAP 1", "UMAP 1",
      "UMAP1", "UMAP1", "HARMONY 1", "HARMONY 1", "HARMONY 1", "HARMONY 1"
    )
  )
})

test_that("`compact_num` works as expected", {
  expect_equal(
    compact_num(c(1, 12, 123, 1234, 12345, 123456, 1234567, 12345678, 123456789)),
    c(
      "1.0", "12.0", "123.0", "1.2k", "12.3k", "123.5k", "1.2M",
      "12.3M", "123.5M"
    )
  )
})
