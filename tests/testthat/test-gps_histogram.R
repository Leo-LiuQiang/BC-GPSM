test_that("gps_histogram returns a ggplot object for valid GPS columns", {
  dat <- data.frame(
    gps_A = c(0.20, 0.35, 0.40, 0.25, 0.30),
    gps_B = c(0.50, 0.40, 0.30, 0.45, 0.35),
    gps_C = c(0.30, 0.25, 0.30, 0.30, 0.35)
  )

  plot <- gps_histogram(dat, bins = 5)

  expect_s3_class(plot, "ggplot")
  expect_equal(plot$labels$x, "logit(gps)")
})

test_that("gps_histogram requires at least two GPS columns", {
  expect_error(gps_histogram(data.frame(gps_A = c(0.2, 0.3))), "two or more")
})
