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

test_that("`order_sample_alias_factors` works as expected", {
  a_tibble <-
    tibble(
      sample_alias = letters,
      another_alias = letters
    ) %>%
    mutate(another_column = 1:26)


  expect_equal(
    order_sample_alias_factors(a_tibble, levels = rev(letters)),
    structure(list(
      sample_alias = structure(26:1, levels = c(
        "z",
        "y", "x", "w", "v", "u", "t", "s", "r", "q", "p", "o", "n", "m",
        "l", "k", "j", "i", "h", "g", "f", "e", "d", "c", "b", "a"
      ), class = "factor"),
      another_alias = c(
        "a", "b", "c", "d", "e", "f", "g", "h",
        "i", "j", "k", "l", "m", "n", "o", "p", "q", "r", "s", "t",
        "u", "v", "w", "x", "y", "z"
      ), another_column = 1:26
    ), row.names = c(
      NA,
      -26L
    ), class = c("tbl_df", "tbl", "data.frame"))
  )

  expect_no_error(order_sample_alias_factors(a_tibble, levels = rev(letters), column_name = "another_alias"))
  expect_error(order_sample_alias_factors(a_tibble, levels = letters[1:4]))
  expect_error(order_sample_alias_factors(a_tibble, levels = LETTERS))
})
