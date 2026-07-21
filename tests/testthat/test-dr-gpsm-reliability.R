reliability_args <- function(data, ...) {
  utils::modifyList(
    list(
        data = data,
        treatment = 1,
        covariate = 2:3,
        outcome = 4,
        gps_model = "logit",
        outcome_model = "lm",
        folds = 2,
        fold_seed = 3101,
        gps_seed = 3102,
        outcome_seed = 3103,
        nboot = 20
      ),
    list(...),
    keep.null = TRUE
  )
}

test_that("dr_gpsm default reference equals an explicit last treatment level", {
  dat <- make_three_arm_data(n_per_arm = 24, seed = 301)

  default_fit <- do.call(dr_gpsm, reliability_args(dat))
  explicit_fit <- do.call(
    dr_gpsm,
    reliability_args(dat, treatment_ref = "C")
  )

  expect_named(default_fit$estimate, c("AvC", "BvC", "BvA"))
  expect_equal(default_fit, explicit_fit, tolerance = 1e-12)
})

test_that("dr_gpsm honors an explicit reference treatment", {
  dat <- make_strong_effect_data(n_per_arm = 30, seed = 302)
  fit <- do.call(dr_gpsm, reliability_args(dat, treatment_ref = "A"))

  expect_named(fit$estimate, c("BvA", "CvA", "CvB"))
  expect_true(all(fit$estimate > 0))
})

test_that("dr_gpsm supports numeric binary outcomes", {
  dat <- make_binary_three_arm_data(n_per_arm = 30, seed = 303)
  fit <- do.call(dr_gpsm, reliability_args(dat, treatment_ref = "A"))

  expect_true(all(is.finite(fit$estimate)))
  expect_true(all(is.finite(fit$ci_lower)))
  expect_true(all(is.finite(fit$ci_upper)))
  expect_true(all(fit$ci_lower <= fit$ci_upper))
})

test_that("dr_gpsm rejects missing outcomes with an actionable error", {
  dat <- make_three_arm_data(n_per_arm = 20, seed = 304)
  dat$y[c(2, 19, 41)] <- NA_real_

  expect_error(
    do.call(dr_gpsm, reliability_args(dat)),
    "outcome.*3 missing values.*remove or impute",
    ignore.case = TRUE
  )
})

test_that("dr_gpsm runs two-step calibration deterministically", {
  skip_if_not_installed("glmnet")
  dat <- make_three_arm_data(n_per_arm = 30, seed = 305)
  args <- reliability_args(
    dat,
    treatment_ref = "A",
    two_step_calibration = TRUE,
    outcome_params = list(two_step = list(lambda = 0.1))
  )

  fit_1 <- do.call(dr_gpsm, args)
  fit_2 <- do.call(dr_gpsm, args)

  expect_equal(fit_1, fit_2, tolerance = 1e-12)
  expect_true(all(is.finite(fit_1$estimate)))
})

test_that("dr_gpsm supports match_ratio greater than one", {
  dat <- make_three_arm_data(n_per_arm = 24, seed = 306)
  fit <- do.call(dr_gpsm, reliability_args(dat, match_ratio = 2L))

  expect_true(all(is.finite(fit$estimate)))
  expect_true(all(fit$ci_lower <= fit$ci_upper))
})

test_that("dr_gpsm is reproducible under fixed seeds", {
  dat <- make_three_arm_data(n_per_arm = 24, seed = 307)
  args <- reliability_args(dat, treatment_ref = "A")

  fit_1 <- do.call(dr_gpsm, args)
  fit_2 <- do.call(dr_gpsm, args)

  expect_equal(fit_1, fit_2, tolerance = 1e-12)
})

test_that("dr_gpsm supports factor covariates", {
  dat <- make_three_arm_data(n_per_arm = 24, seed = 308)
  dat$site <- factor(rep(c("North", "South", "West"), length.out = nrow(dat)))
  dat <- dat[, c("trt", "x1", "x2", "site", "y")]

  args <- reliability_args(dat)
  args$covariate <- 2:4
  args$outcome <- 5
  fit <- do.call(dr_gpsm, args)

  expect_true(all(is.finite(fit$estimate)))
})

test_that("dr_gpsm handles viable rare arms and rejects impossible folds", {
  viable <- make_rare_arm_data(n_common = 24, n_rare = 8, seed = 309)
  fit <- do.call(dr_gpsm, reliability_args(viable))
  expect_true(all(is.finite(fit$estimate)))

  too_rare <- make_rare_arm_data(n_common = 12, n_rare = 1, seed = 310)
  expect_error(
    do.call(dr_gpsm, reliability_args(too_rare)),
    "treatment arm.*C.*1 observation.*folds = 2.*reduce.*folds",
    ignore.case = TRUE
  )
})

test_that("bootstrap intervals are ordered and preserve expected contrast signs", {
  dat <- make_strong_effect_data(n_per_arm = 35, seed = 311)
  fit <- do.call(
    dr_gpsm,
    reliability_args(dat, treatment_ref = "A", nboot = 50)
  )

  expect_named(fit$estimate, c("BvA", "CvA", "CvB"))
  expect_true(all(fit$estimate > 0))
  expect_true(all(is.finite(fit$ci_lower)))
  expect_true(all(is.finite(fit$ci_upper)))
  expect_true(all(fit$ci_lower < fit$ci_upper))
})
