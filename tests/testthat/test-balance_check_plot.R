test_that("balance_check_plot returns a plot and balance table for covariate matching", {
  dat <- make_three_arm_data(n_per_arm = 8, seed = 41)

  out <- balance_check_plot(
    data = dat,
    treatment = 1,
    covariate = 2:3,
    match_on = "covariates",
    return_data = TRUE
  )

  expect_s3_class(out$plot, "ggplot")
  expect_s3_class(out$balance, "data.frame")
  expect_s3_class(out$plot_data, "data.frame")
  expect_named(out$balance, c("contrast", "covariate", "sample", "abs_smd"))
  expect_true(all(out$balance$sample %in% c("Before matching", "After matching")))
  expect_true(all(as.character(out$plot_data$sample) %in% c("All", "Matched")))
  expect_true(all(is.finite(out$balance$abs_smd)))
})

test_that("balance_check_plot uses supplied loggps columns", {
  dat <- make_three_arm_data(n_per_arm = 8, seed = 42)
  dat <- gps_pre_process(
    dat,
    treatment = 1,
    treatment_ref = "A",
    covariate = 2:3,
    fit_gps = FALSE
  )
  dat$loggps_B <- seq(-1, 1, length.out = nrow(dat))
  dat$loggps_C <- seq(1, -1, length.out = nrow(dat))

  out <- balance_check_plot(
    data = dat,
    treatment = 1,
    treatment_ref = "A",
    covariate = 2:3,
    match_on = "gps",
    fit_gps = FALSE,
    style = "faceted",
    return_data = TRUE
  )

  expect_s3_class(out$plot, "ggplot")
  expect_s3_class(out$plot_data, "data.frame")
  expect_equal(sort(unique(as.character(out$balance$sample))),
               c("After matching", "Before matching"))
})

test_that("balance_check_plot validates missing loggps columns", {
  dat <- make_three_arm_data(n_per_arm = 8, seed = 43)

  expect_error(
    balance_check_plot(
      data = dat,
      treatment = 1,
      treatment_ref = "A",
      covariate = 2:3,
      match_on = "gps",
      fit_gps = FALSE
    ),
    "Missing log-GPS columns"
  )
})
