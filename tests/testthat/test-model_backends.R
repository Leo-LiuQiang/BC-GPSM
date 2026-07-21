test_that("new GPS model choices can be used without fitting diagnostics", {
  dat <- make_three_arm_data(n_per_arm = 30)

  xgb <- gps_pre_process(
    data = dat,
    treatment = 1,
    treatment_ref = "A",
    covariate = 2:3,
    gps_model = "xgboost",
    fit_gps = FALSE
  )
  rng <- gps_pre_process(
    data = dat,
    treatment = 1,
    treatment_ref = "A",
    covariate = 2:3,
    gps_model = "ranger",
    fit_gps = FALSE
  )

  expect_s3_class(xgb, "data.frame")
  expect_s3_class(rng, "data.frame")
  expect_true(all(grepl("^gps_", names(xgb)[grep("^gps_", names(xgb))])))
  expect_true(all(grepl("^gps_", names(rng)[grep("^gps_", names(rng))])))
})

test_that("xgboost GPS backend either runs or reports missing dependency clearly", {
  dat <- make_three_arm_data(n_per_arm = 30)
  train_id <- c(1:20, 31:50, 61:80)

  call_backend <- function() {
    pmatch_fit_ps_fold(
      train_df = dat[train_id, ],
      eval_df = dat[-train_id, ],
      t_name = "trt",
      x_names = c("x1", "x2"),
      model = "xgboost",
      tr_levels = levels(dat$trt),
      gps_params = list(xgboost = list(nrounds = 2L, nthread = 1L, verbose = 0))
    )
  }

  if (!requireNamespace("xgboost", quietly = TRUE)) {
    expect_error(call_backend(), "Package 'xgboost' is required")
  } else {
    out <- call_backend()
    expect_equal(dim(out$e_hat_eval), c(nrow(dat) - length(train_id), nlevels(dat$trt)))
    expect_equal(rowSums(out$e_hat_eval), rep(1, nrow(out$e_hat_eval)), tolerance = 1e-6)
  }
})

test_that("ranger GPS backend either runs or reports missing dependency clearly", {
  dat <- make_three_arm_data(n_per_arm = 30)
  train_id <- c(1:20, 31:50, 61:80)

  call_backend <- function() {
    pmatch_fit_ps_fold(
      train_df = dat[train_id, ],
      eval_df = dat[-train_id, ],
      t_name = "trt",
      x_names = c("x1", "x2"),
      model = "ranger",
      tr_levels = levels(dat$trt),
      gps_params = list(ranger = list(num.trees = 5L, num.threads = 1L))
    )
  }

  if (!requireNamespace("ranger", quietly = TRUE)) {
    expect_error(call_backend(), "Package 'ranger' is required")
  } else {
    out <- call_backend()
    expect_equal(dim(out$e_hat_eval), c(nrow(dat) - length(train_id), nlevels(dat$trt)))
    expect_equal(rowSums(out$e_hat_eval), rep(1, nrow(out$e_hat_eval)), tolerance = 1e-6)
  }
})
