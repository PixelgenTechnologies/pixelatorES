test_that("create_sample_palette works as expected", {
  expect_no_error(create_sample_palette("letters"))
  expect_no_error(create_sample_palette(letters[1:2]))
  expect_no_error(create_sample_palette(rep(letters[1:2], 4)))
  expect_no_error(create_sample_palette(letters[1:6]))
  expect_no_error(create_sample_palette(letters))
})
