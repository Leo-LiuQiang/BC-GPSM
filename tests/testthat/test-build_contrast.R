test_that("build_contrast uses the last level as the default reference", {
  contrast <- build_contrast(c("A", "B", "C"))

  expect_equal(dim(contrast), c(3L, 3L))
  expect_equal(colnames(contrast), c("A", "B", "C"))
  expect_true(all(c("BvA", "AvC", "BvC") %in% rownames(contrast)))
  expect_equal(unname(contrast["AvC", ]), c(1, 0, -1))
  expect_equal(unname(contrast["BvC", ]), c(0, 1, -1))
})

test_that("build_contrast orients contrasts against an explicit reference", {
  contrast <- build_contrast(c("A", "B", "C"), ref = "A")

  expect_equal(unname(contrast["BvA", ]), c(-1, 1, 0))
  expect_equal(unname(contrast["CvA", ]), c(-1, 0, 1))
  expect_equal(unname(contrast["CvB", ]), c(0, -1, 1))
})

test_that("build_contrast validates treatment levels and reference", {
  expect_error(build_contrast("A"), "at least two")
  expect_error(build_contrast(c("A", "B"), ref = "C"), "ref")
})
