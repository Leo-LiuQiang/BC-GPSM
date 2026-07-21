test_that("pmatch_compute_V uses the last probability column by default", {
  e_mat <- matrix(
    c(0.20, 0.50, 0.30,
      0.40, 0.40, 0.20),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(NULL, c("A", "B", "C"))
  )

  V <- pmatch_compute_V(e_mat)

  expect_equal(dim(V), c(2L, 2L))
  expect_equal(colnames(V), c("A_vs_C", "B_vs_C"))
  expect_equal(unname(V[1, ]), log(c(0.20, 0.50) / 0.30))
})

test_that("pmatch_compute_V supports reference column names and indices", {
  e_mat <- matrix(
    c(0.20, 0.50, 0.30,
      0.40, 0.40, 0.20),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(NULL, c("A", "B", "C"))
  )

  by_name <- pmatch_compute_V(e_mat, ref_col = "A")
  by_index <- pmatch_compute_V(e_mat, ref_col = 1)

  expect_equal(by_name, by_index)
  expect_equal(colnames(by_name), c("B_vs_A", "C_vs_A"))
  expect_equal(unname(by_name[2, ]), log(c(0.40, 0.20) / 0.40))
})

test_that("pmatch_compute_V validates reference columns", {
  e_mat <- matrix(c(0.4, 0.6, 0.5, 0.5), ncol = 2, byrow = TRUE)

  expect_error(pmatch_compute_V(e_mat, ref_col = 3), "ref_col")
  expect_error(pmatch_compute_V(e_mat, ref_col = "missing"), "ref_col")
})
