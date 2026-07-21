test_that("dr_gpsm reports invalid data and column selections clearly", {
  dat <- make_three_arm_data(n_per_arm = 8, seed = 401)

  expect_error(
    dr_gpsm(as.matrix(dat), 1, covariate = 2:3, outcome = 4),
    "data.*data.frame"
  )
  expect_error(
    dr_gpsm(dat, 0, covariate = 2:3, outcome = 4),
    "treatment.*outside 1:4"
  )
  expect_error(
    dr_gpsm(dat, 1, covariate = c(2, 2), outcome = 4),
    "covariate.*duplicate"
  )
  expect_error(
    dr_gpsm(dat, 1, covariate = 2:4, outcome = 4),
    "distinct columns"
  )
})

test_that("dr_gpsm reports invalid computation settings clearly", {
  dat <- make_three_arm_data(n_per_arm = 8, seed = 402)

  expect_error(
    dr_gpsm(dat, 1, covariate = 2:3, outcome = 4, folds = 1),
    "folds.*greater than or equal to 2"
  )
  expect_error(
    dr_gpsm(dat, 1, covariate = 2:3, outcome = 4, nboot = 0),
    "nboot.*greater than or equal to 1"
  )
  expect_error(
    dr_gpsm(dat, 1, covariate = 2:3, outcome = 4, match_ratio = 1.5),
    "match_ratio.*integer"
  )
  expect_error(
    dr_gpsm(
      dat, 1, covariate = 2:3, outcome = 4,
      outcome_model = "none", two_step_calibration = TRUE
    ),
    "requires an outcome model"
  )
})

test_that("dr_gpsm explains unsupported outcome encodings", {
  dat <- make_binary_three_arm_data(n_per_arm = 10, seed = 403)
  dat$y <- factor(dat$y, levels = c(0, 1), labels = c("No", "Yes"))

  expect_error(
    dr_gpsm(dat, 1, covariate = 2:3, outcome = 4, outcome_model = "lm"),
    "requires a numeric outcome.*0/1"
  )
})

test_that("dr_gpsm checks matching donors before model fitting", {
  dat <- make_rare_arm_data(n_common = 12, n_rare = 2, seed = 404)

  expect_error(
    dr_gpsm(
      dat, 1, covariate = 2:3, outcome = 4,
      folds = 2, match_ratio = 3, nboot = 2
    ),
    "arm 'C'.*fewer than.*match_ratio = 3.*Reduce"
  )
})

test_that("dr_gpsm reports invalid reference treatments clearly", {
  dat <- make_three_arm_data(n_per_arm = 8, seed = 405)

  expect_error(
    dr_gpsm(
      dat, 1, treatment_ref = "Control", covariate = 2:3,
      outcome = 4, nboot = 2
    ),
    "treatment_ref.*level"
  )
})
