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

make_strong_effect_data <- function(n_per_arm = 30, seed = 1) {
  set.seed(seed)
  trt <- factor(rep(c("A", "B", "C"), each = n_per_arm),
                levels = c("A", "B", "C"))
  n <- length(trt)
  x1 <- rnorm(n)
  x2 <- runif(n, -1, 1)
  y <- 0.8 * (trt == "B") +
    1.8 * (trt == "C") +
    0.2 * x1 -
    0.1 * x2 +
    rnorm(n, sd = 0.08)

  data.frame(trt = trt, x1 = x1, x2 = x2, y = y)
}

make_binary_three_arm_data <- function(n_per_arm = 30, seed = 1) {
  dat <- make_three_arm_data(n_per_arm = n_per_arm, seed = seed)
  eta <- -0.8 +
    0.7 * (dat$trt == "B") +
    1.4 * (dat$trt == "C") +
    0.25 * dat$x1 -
    0.15 * dat$x2

  set.seed(seed + 1000L)
  dat$y <- stats::rbinom(nrow(dat), size = 1L, prob = stats::plogis(eta))
  dat
}

make_rare_arm_data <- function(n_common = 24, n_rare = 8, seed = 1) {
  stopifnot(n_common >= 1L, n_rare >= 1L)
  set.seed(seed)
  trt <- factor(
    c(rep("A", n_common), rep("B", n_common), rep("C", n_rare)),
    levels = c("A", "B", "C")
  )
  n <- length(trt)
  x1 <- rnorm(n)
  x2 <- runif(n)
  y <- 0.5 * (trt == "B") +
    1.0 * (trt == "C") +
    0.3 * x1 -
    0.2 * x2 +
    rnorm(n, sd = 0.2)

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
