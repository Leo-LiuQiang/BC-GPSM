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
