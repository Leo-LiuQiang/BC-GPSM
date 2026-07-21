test_that("gps_matching returns estimates and per-unit contributions for GPS matching", {
  dat <- make_three_arm_data(n_per_arm = 5, seed = 31)
  dat <- gps_pre_process(
    dat,
    treatment = 1,
    treatment_ref = "A",
    covariate = 2:3,
    fit_gps = FALSE
  )
  dat$loggps_B <- seq(-1, 1, length.out = nrow(dat))
  dat$loggps_C <- seq(1, -1, length.out = nrow(dat))

  contrast <- build_contrast(levels(dat$trt), ref = "A")
  pred <- make_prediction_matrix(dat)

  fit <- gps_matching(
    data = dat,
    treatment = 1,
    outcome = 4,
    pred = pred,
    contrast = contrast,
    return_tau = TRUE,
    do_boot = FALSE
  )

  expect_named(fit$estimate, rownames(contrast))
  expect_true(all(is.finite(fit$estimate)))
  expect_true(all(is.na(fit$ci_lower)))
  expect_true(all(is.na(fit$ci_upper)))
  expect_equal(dim(fit$tau), c(nrow(contrast), nrow(dat)))
  expect_equal(dim(fit$psi), c(length(levels(dat$trt)), nrow(dat)))
})

test_that("gps_matching supports covariate matching", {
  dat <- make_three_arm_data(n_per_arm = 5, seed = 32)
  dat <- gps_pre_process(
    dat,
    treatment = 1,
    treatment_ref = "A",
    covariate = 2:3,
    fit_gps = FALSE
  )

  contrast <- build_contrast(levels(dat$trt), ref = "A")
  pred <- make_prediction_matrix(dat)

  fit <- gps_matching(
    data = dat,
    treatment = 1,
    outcome = 4,
    pred = pred,
    contrast = contrast,
    match_on = "covariates",
    covariate = 2:3,
    do_boot = FALSE
  )

  expect_named(fit$estimate, rownames(contrast))
  expect_true(all(is.finite(fit$estimate)))
})

test_that("gps_matching validates required inputs", {
  dat <- make_three_arm_data(n_per_arm = 4, seed = 33)
  dat <- gps_pre_process(
    dat,
    treatment = 1,
    treatment_ref = "A",
    covariate = 2:3,
    fit_gps = FALSE
  )
  contrast <- build_contrast(levels(dat$trt), ref = "A")
  pred <- make_prediction_matrix(dat)

  expect_error(
    gps_matching(
      data = dat,
      treatment = 1,
      outcome = 4,
      pred = pred,
      contrast = contrast,
      match_on = "covariates",
      do_boot = FALSE
    ),
    "covariate"
  )
})
