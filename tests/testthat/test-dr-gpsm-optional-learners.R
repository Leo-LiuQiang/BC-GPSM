optional_fit_args <- function(data, ...) {
  utils::modifyList(
    list(
      data = data,
      treatment = 1,
      treatment_ref = "A",
      covariate = 2:3,
      outcome = 4,
      gps_model = "logit",
      outcome_model = "lm",
      folds = 2,
      fold_seed = 5101,
      gps_seed = 5102,
      outcome_seed = 5103,
      nboot = 3
    ),
    list(...),
    keep.null = TRUE
  )
}

expect_valid_optional_fit <- function(fit) {
  expect_named(fit$estimate, c("BvA", "CvA", "CvB"))
  expect_true(all(is.finite(fit$estimate)))
  expect_true(all(is.finite(fit$ci_lower)))
  expect_true(all(is.finite(fit$ci_upper)))
}

test_that("dr_gpsm completes a full GBM GPS run when installed", {
  skip_if_not_installed("gbm")
  dat <- make_three_arm_data(n_per_arm = 24, seed = 501)
  fit <- do.call(
    dr_gpsm,
    optional_fit_args(
      dat,
      gps_model = "gbm",
      gps_params = list(gbm = list(
        n.trees = 20L,
        interaction.depth = 1L,
        shrinkage = 0.1,
        bag.fraction = 0.8,
        cv.folds = 2L,
        n.minobsinnode = 5L,
        n.cores = 1L,
        verbose = FALSE
      ))
    )
  )
  expect_valid_optional_fit(fit)
})

test_that("dr_gpsm completes a full GAM GPS run when installed", {
  skip_if_not_installed("mgcv")
  dat <- make_three_arm_data(n_per_arm = 24, seed = 502)
  fit <- do.call(
    dr_gpsm,
    optional_fit_args(
      dat,
      gps_model = "gam",
      gps_params = list(gam = list(df_max = 3L, method = "REML"))
    )
  )
  expect_valid_optional_fit(fit)
})

test_that("dr_gpsm completes a full XGBoost GPS run when installed", {
  skip_if_not_installed("xgboost")
  dat <- make_three_arm_data(n_per_arm = 24, seed = 503)
  fit <- do.call(
    dr_gpsm,
    optional_fit_args(
      dat,
      gps_model = "xgboost",
      gps_params = list(xgboost = list(
        nrounds = 3L,
        max_depth = 1L,
        eta = 0.2,
        subsample = 1,
        colsample_bytree = 1,
        nthread = 1L,
        verbose = 0
      ))
    )
  )
  expect_valid_optional_fit(fit)
})

test_that("dr_gpsm completes a full ranger GPS run when installed", {
  skip_if_not_installed("ranger")
  dat <- make_three_arm_data(n_per_arm = 24, seed = 504)
  fit <- do.call(
    dr_gpsm,
    optional_fit_args(
      dat,
      gps_model = "ranger",
      gps_params = list(ranger = list(
        num.trees = 10L,
        min.node.size = 3L,
        sample.fraction = 1,
        replace = TRUE,
        num.threads = 1L
      ))
    )
  )
  expect_valid_optional_fit(fit)
})

test_that("dr_gpsm completes a full random forest outcome run when installed", {
  skip_if_not_installed("randomForest")
  dat <- make_three_arm_data(n_per_arm = 24, seed = 505)
  fit <- do.call(
    dr_gpsm,
    optional_fit_args(
      dat,
      outcome_model = "rf",
      outcome_params = list(rf = list(ntree = 10L, mtry = 1L, nodesize = 3L))
    )
  )
  expect_valid_optional_fit(fit)
})

test_that("dr_gpsm completes a full GBM outcome run when installed", {
  skip_if_not_installed("gbm")
  dat <- make_three_arm_data(n_per_arm = 30, seed = 506)
  fit <- do.call(
    dr_gpsm,
    optional_fit_args(
      dat,
      outcome_model = "gbm",
      outcome_params = list(gbm = list(
        n.trees = 20L,
        interaction.depth = 1L,
        shrinkage = 0.1,
        bag.fraction = 0.8,
        cv.folds = 2L,
        n.minobsinnode = 2L,
        n.cores = 1L,
        verbose = FALSE
      ))
    )
  )
  expect_valid_optional_fit(fit)
})

test_that("dr_gpsm completes a full GAM outcome run when installed", {
  skip_if_not_installed("mgcv")
  dat <- make_three_arm_data(n_per_arm = 24, seed = 507)
  fit <- do.call(
    dr_gpsm,
    optional_fit_args(
      dat,
      outcome_model = "gam",
      outcome_params = list(gam = list(df_max = 3L, method = "REML"))
    )
  )
  expect_valid_optional_fit(fit)
})

test_that("dr_gpsm completes a full XGBoost outcome run when installed", {
  skip_if_not_installed("xgboost")
  dat <- make_three_arm_data(n_per_arm = 24, seed = 508)
  fit <- do.call(
    dr_gpsm,
    optional_fit_args(
      dat,
      outcome_model = "xgboost",
      outcome_params = list(xgboost = list(
        nrounds = 3L,
        max_depth = 1L,
        eta = 0.2,
        subsample = 1,
        colsample_bytree = 1,
        nthread = 1L,
        verbose = 0
      ))
    )
  )
  expect_valid_optional_fit(fit)
})

test_that("dr_gpsm completes a full ranger outcome run when installed", {
  skip_if_not_installed("ranger")
  dat <- make_three_arm_data(n_per_arm = 24, seed = 509)
  fit <- do.call(
    dr_gpsm,
    optional_fit_args(
      dat,
      outcome_model = "ranger",
      outcome_params = list(ranger = list(
        num.trees = 10L,
        min.node.size = 3L,
        sample.fraction = 1,
        replace = TRUE,
        num.threads = 1L
      ))
    )
  )
  expect_valid_optional_fit(fit)
})
