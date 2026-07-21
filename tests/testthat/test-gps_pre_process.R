test_that("gps_pre_process creates placeholders and defaults to the last reference level", {
  dat <- make_three_arm_data(n_per_arm = 4, seed = 11)
  dat$x1[1] <- NA

  out <- gps_pre_process(
    data = dat,
    treatment = 1,
    covariate = 2:3,
    fit_gps = FALSE
  )

  expect_equal(nrow(out), nrow(dat) - 1L)
  expect_equal(levels(out$trt), c("C", "A", "B"))
  expect_true(all(c("gps_C", "gps_A", "gps_B", "loggps_A", "loggps_B") %in% names(out)))
  expect_true(all(is.na(out$gps_C)))
  expect_true(all(is.na(out$loggps_A)))
})

test_that("gps_pre_process honors explicit treatment_ref", {
  dat <- make_three_arm_data(n_per_arm = 4, seed = 12)

  out <- gps_pre_process(
    data = dat,
    treatment = 1,
    treatment_ref = "B",
    covariate = 2:3,
    fit_gps = FALSE
  )

  expect_equal(levels(out$trt), c("B", "A", "C"))
  expect_true(all(c("gps_B", "gps_A", "gps_C", "loggps_A", "loggps_C") %in% names(out)))
  expect_error(
    gps_pre_process(dat, treatment = 1, treatment_ref = "D", covariate = 2:3),
    "treatment_ref"
  )
})

test_that("gps_pre_process can fit multinomial logit GPS diagnostics", {
  dat <- make_three_arm_data(n_per_arm = 12, seed = 13)

  out <- gps_pre_process(
    data = dat,
    treatment = 1,
    treatment_ref = "A",
    covariate = 2:3,
    gps_model = "logit",
    fit_gps = TRUE
  )

  gps_cols <- grep("^gps_", names(out), value = TRUE)
  loggps_cols <- grep("^loggps_", names(out), value = TRUE)

  expect_false(anyNA(out[gps_cols]))
  expect_equal(unname(rowSums(out[gps_cols])), rep(1, nrow(out)), tolerance = 1e-6)
  expect_true(all(is.finite(as.matrix(out[loggps_cols]))))
})

expect_valid_gps_diagnostics <- function(out) {
  gps_cols <- grep("^gps_", names(out), value = TRUE)
  loggps_cols <- grep("^loggps_", names(out), value = TRUE)

  expect_length(gps_cols, 3L)
  expect_length(loggps_cols, 2L)
  expect_false(anyNA(out[gps_cols]))
  expect_equal(unname(rowSums(out[gps_cols])), rep(1, nrow(out)),
               tolerance = 1e-5)
  expect_true(all(is.finite(as.matrix(out[loggps_cols]))))
}

test_that("gps_pre_process reports invalid inputs and empty analyses clearly", {
  dat <- make_three_arm_data(n_per_arm = 4, seed = 14)

  expect_error(
    gps_pre_process(as.matrix(dat), 1, covariate = 2:3),
    "data.*data.frame"
  )
  expect_error(
    gps_pre_process(dat, treatment = c(1, 2), covariate = 2:3),
    "treatment.*single numeric"
  )
  expect_error(
    gps_pre_process(dat, treatment = 1, covariate = character()),
    "covariate.*numeric vector"
  )
  expect_error(
    gps_pre_process(dat, treatment = 1, covariate = c(2, 8)),
    "indices out of range"
  )

  no_complete <- dat
  no_complete$x1 <- NA_real_
  expect_error(
    gps_pre_process(no_complete, treatment = 1, covariate = 2:3),
    "No complete cases"
  )

  one_arm <- droplevels(dat[dat$trt == "A", ])
  expect_error(
    gps_pre_process(one_arm, treatment = 1, covariate = 2:3),
    "at least 2 levels"
  )
})

test_that("gps_pre_process fits GBM GPS diagnostics when installed", {
  skip_if_not_installed("gbm")
  dat <- make_three_arm_data(n_per_arm = 15, seed = 15)

  out <- gps_pre_process(
    dat, treatment = 1, treatment_ref = "A", covariate = 2:3,
    gps_model = "gbm", fit_gps = TRUE,
    gps_params = list(gbm = list(
      n.trees = 12L, interaction.depth = 1L, shrinkage = 0.1,
      cv.folds = 0L, n.minobsinnode = 3L, n.cores = 1L,
      verbose = FALSE
    ))
  )

  expect_valid_gps_diagnostics(out)
})

test_that("gps_pre_process fits XGBoost GPS diagnostics when installed", {
  skip_if_not_installed("xgboost")
  dat <- make_three_arm_data(n_per_arm = 15, seed = 16)

  out <- gps_pre_process(
    dat, treatment = 1, treatment_ref = "A", covariate = 2:3,
    gps_model = "xgboost", fit_gps = TRUE,
    gps_params = list(xgboost = list(
      nrounds = 3L, max_depth = 1L, eta = 0.2, subsample = 1,
      colsample_bytree = 1, nthread = 1L, verbose = 0
    ))
  )

  expect_valid_gps_diagnostics(out)
})

test_that("gps_pre_process fits ranger GPS diagnostics when installed", {
  skip_if_not_installed("ranger")
  dat <- make_three_arm_data(n_per_arm = 15, seed = 17)

  out <- gps_pre_process(
    dat, treatment = 1, treatment_ref = "A", covariate = 2:3,
    gps_model = "ranger", fit_gps = TRUE,
    gps_params = list(ranger = list(
      num.trees = 12L, min.node.size = 3L, sample.fraction = 1,
      replace = TRUE, num.threads = 1L
    ))
  )

  expect_valid_gps_diagnostics(out)
})

test_that("gps_pre_process fits GAM GPS diagnostics when installed", {
  skip_if_not_installed("VGAM")
  dat <- make_three_arm_data(n_per_arm = 15, seed = 18)
  dat$site <- factor(rep(c("North", "South", "West"), length.out = nrow(dat)))

  out <- gps_pre_process(
    dat, treatment = 1, treatment_ref = "A", covariate = 2:3,
    gps_model = "gam", fit_gps = TRUE
  )

  expect_valid_gps_diagnostics(out)
})

test_that("gps_pre_process explains unsupported automatic tuning", {
  dat <- make_three_arm_data(n_per_arm = 5, seed = 19)

  expect_error(
    gps_pre_process(
      dat, treatment = 1, covariate = 2:3, gps_model = "xgboost",
      tune = TRUE
    ),
    "not currently implemented.*gps_params\\$xgboost"
  )
  expect_error(
    gps_pre_process(
      dat, treatment = 1, covariate = 2:3, gps_model = "ranger",
      tune = TRUE
    ),
    "not currently implemented.*gps_params\\$ranger"
  )
})
