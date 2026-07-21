make_three_arm_data <- function(n_per_arm = 20, seed = 1) {
  set.seed(seed)
  trt <- factor(rep(c("A", "B", "C"), each = n_per_arm),
                levels = c("A", "B", "C"))
  n <- length(trt)
  x1 <- rnorm(n)
  x2 <- runif(n)
  y <- 1 +
    0.4 * (trt == "B") +
    0.8 * (trt == "C") +
    0.6 * x1 -
    0.3 * x2 +
    rnorm(n, sd = 0.4)

  data.frame(trt = trt, x1 = x1, x2 = x2, y = y)
}

make_prediction_matrix <- function(data, treatment_col = "trt", outcome_col = "y") {
  trt_levels <- levels(data[[treatment_col]])
  arm_means <- tapply(data[[outcome_col]], data[[treatment_col]], mean)[trt_levels]

  matrix(
    rep(as.numeric(arm_means), each = nrow(data)),
    nrow = length(trt_levels),
    byrow = TRUE,
    dimnames = list(trt_levels, NULL)
  )
}
