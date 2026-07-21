test_that("dr_gpsm runs a small cross-fitted workflow", {
  dat <- make_three_arm_data(n_per_arm = 20, seed = 41)

  fit <- dr_gpsm(
    data = dat,
    treatment = 1,
    treatment_ref = "A",
    covariate = 2:3,
    outcome = 4,
    gps_model = "logit",
    outcome_model = "lm",
    folds = 2,
    fold_seed = 42,
    gps_seed = 43,
    outcome_seed = 44,
    nboot = 5
  )

  expect_equal(names(fit), c("estimate", "ci_lower", "ci_upper"))
  expect_named(fit$estimate, c("BvA", "CvA", "CvB"))
  expect_true(all(is.finite(fit$estimate)))
  expect_true(all(is.finite(fit$ci_lower)))
  expect_true(all(is.finite(fit$ci_upper)))
  expect_true(all(fit$ci_lower <= fit$ci_upper))
})

test_that("dr_gpsm requires hist_path when saving a GPS histogram", {
  dat <- make_three_arm_data(n_per_arm = 6, seed = 45)

  expect_error(
    dr_gpsm(
      data = dat,
      treatment = 1,
      treatment_ref = "A",
      covariate = 2:3,
      outcome = 4,
      gps_model = "logit",
      outcome_model = "lm",
      folds = 2,
      nboot = 2,
      hist = TRUE
    ),
    "hist_path"
  )
})
