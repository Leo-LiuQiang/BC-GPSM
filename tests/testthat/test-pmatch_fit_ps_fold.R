test_that("pmatch_fit_ps_fold returns probabilities in the requested treatment order", {
  dat <- make_three_arm_data(n_per_arm = 15, seed = 21)
  ids_by_arm <- split(seq_len(nrow(dat)), dat$trt)
  train_id <- unlist(lapply(ids_by_arm, head, 10), use.names = FALSE)
  eval_id <- setdiff(seq_len(nrow(dat)), train_id)

  fit <- pmatch_fit_ps_fold(
    train_df = dat[train_id, ],
    eval_df = dat[eval_id, ],
    t_name = "trt",
    x_names = c("x1", "x2"),
    model = "logit",
    tr_levels = c("C", "A", "B"),
    seed = 22
  )

  expect_equal(fit$levels, c("C", "A", "B"))
  expect_equal(colnames(fit$e_hat_eval), c("C", "A", "B"))
  expect_equal(nrow(fit$e_hat_eval), length(eval_id))
  expect_equal(unname(rowSums(fit$e_hat_eval)), rep(1, length(eval_id)), tolerance = 1e-6)
  expect_true(all(is.finite(fit$e_hat_eval)))
})
