#' Covariate balance check plot
#'
#' Computes pairwise absolute standardized mean differences (SMDs) before and
#' after nearest-neighbor matching, then returns a Love-plot style balance plot.
#' Matching can be performed on existing \code{loggps_*} columns, on a
#' full-sample diagnostic GPS fit, or directly on covariates.
#'
#' @param data Data frame containing treatment and covariates.
#' @param treatment Single numeric column index for treatment.
#' @param covariate Numeric vector of covariate column indices to check.
#' @param treatment_ref Optional reference treatment level. If \code{NULL}, the
#'   last observed treatment level is used when a GPS fit is needed.
#' @param match_on Matching scale: \code{"gps"} or \code{"covariates"}.
#' @param gps_model GPS model used when \code{match_on = "gps"} and
#'   \code{fit_gps = TRUE}; one of \code{"logit"}, \code{"gbm"},
#'   \code{"gam"}, \code{"xgboost"}, or \code{"ranger"}. Flexible learners
#'   require optional packages; XGBoost and ranger are experimental.
#' @param gps_params Optional named list of model-specific GPS parameters.
#' @param fit_gps Logical. If \code{TRUE} and \code{match_on = "gps"}, fit a
#'   full-sample diagnostic GPS model before computing balance. If \code{FALSE},
#'   \code{data} must already contain finite \code{loggps_*} columns.
#' @param cov_distance For covariate matching: \code{"mahalanobis"} or
#'   \code{"euclidean"}.
#' @param standardize If \code{TRUE}, center and scale covariate matching
#'   features before matching.
#' @param ridge Small ridge used for Mahalanobis whitening.
#' @param match_ratio Integer >= 1; number of nearest neighbors per target arm.
#' @param threshold Reference line for acceptable absolute SMD.
#' @param style Plot style. \code{"love"} aggregates across pairwise treatment
#'   contrasts by taking the maximum absolute SMD for each covariate and sample;
#'   \code{"faceted"} shows each pairwise treatment contrast separately.
#' @param return_data If \code{TRUE}, return a list with \code{plot} and
#'   \code{balance}; otherwise return only the plot.
#'
#' @return A \pkg{ggplot2} object, or a list with the plot and balance table
#'   when \code{return_data = TRUE}.
#' @export
#'
#' @examples
#' set.seed(1)
#' n <- 75
#' z1 <- rnorm(n)
#' z2 <- runif(n, -1, 1)
#' dat <- data.frame(
#'   trt = factor(rep(c("A", "B", "C"), each = 25)),
#'   age = 50 + 10 * z1,
#'   biomarker = z2,
#'   severity = 0.5 * z1 - 0.3 * z2 + rnorm(n, sd = 0.7),
#'   prior_tx = rbinom(n, 1, plogis(0.4 * z1)),
#'   site = factor(sample(c("S1", "S2"), n, replace = TRUE))
#' )
#'
#' p <- balance_check_plot(
#'   data = dat,
#'   treatment = 1,
#'   covariate = 2:6,
#'   match_on = "covariates",
#'   style = "love"
#' )
#'
#' p
balance_check_plot <- function(data,
                               treatment,
                               covariate,
                               treatment_ref = NULL,
                               match_on = c("gps", "covariates"),
                               gps_model = c("logit", "gbm", "gam", "xgboost", "ranger"),
                               gps_params = NULL,
                               fit_gps = TRUE,
                               cov_distance = c("mahalanobis", "euclidean"),
                               standardize = FALSE,
                               ridge = 1e-8,
                               match_ratio = 1L,
                               threshold = 0.1,
                               style = c("love", "faceted"),
                               return_data = FALSE) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for balance_check_plot().")
  }
  if (!requireNamespace("RANN", quietly = TRUE)) {
    stop("Package 'RANN' is required for balance_check_plot() matching.")
  }

  match_on <- match.arg(match_on)
  gps_model <- match.arg(gps_model)
  cov_distance <- match.arg(cov_distance)
  style <- match.arg(style)

  if (!is.data.frame(data)) stop("`data` must be a data.frame.")
  if (!is.numeric(treatment) || length(treatment) != 1L) {
    stop("`treatment` must be a single numeric column index.")
  }
  if (!is.numeric(covariate) || length(covariate) < 1L) {
    stop("`covariate` must be a numeric vector of column indices.")
  }
  if (!is.numeric(match_ratio) || length(match_ratio) != 1L || match_ratio < 1) {
    stop("`match_ratio` must be a single integer >= 1.")
  }
  match_ratio <- as.integer(match_ratio)
  if (!is.numeric(threshold) || length(threshold) != 1L || !is.finite(threshold)) {
    stop("`threshold` must be a single finite number.")
  }

  p <- ncol(data)
  if (any(c(treatment, covariate) < 1 | c(treatment, covariate) > p)) {
    stop("treatment/covariate indices out of range.")
  }

  trt_var <- names(data)[treatment]
  cov_vars <- names(data)[covariate]
  cc_cols <- unique(c(trt_var, cov_vars))
  data <- data[stats::complete.cases(data[, cc_cols, drop = FALSE]), , drop = FALSE]
  if (nrow(data) == 0L) stop("No complete cases after filtering on treatment + covariates.")

  if (match_on == "gps" && isTRUE(fit_gps)) {
    data <- gps_pre_process(
      data = data,
      treatment = treatment,
      treatment_ref = treatment_ref,
      covariate = covariate,
      gps_model = gps_model,
      gps_params = gps_params,
      fit_gps = TRUE
    )
  } else {
    trt_fac <- factor(data[[trt_var]])
    original_levels <- levels(trt_fac)
    if (is.null(treatment_ref)) {
      loggps_levels <- sub("^loggps_", "", grep("^loggps_", names(data), value = TRUE))
      inferred_ref <- setdiff(original_levels, loggps_levels)
      treatment_ref <- if (match_on == "gps" && length(inferred_ref) == 1L) {
        inferred_ref
      } else {
        original_levels[length(original_levels)]
      }
    }
    if (!treatment_ref %in% original_levels) {
      stop("`treatment_ref` must be a level in treatment column.")
    }
    data[[trt_var]] <- stats::relevel(trt_fac, ref = treatment_ref)
  }

  trt_fac <- data[[trt_var]]
  if (!is.factor(trt_fac)) stop("Internal error: treatment column is not a factor.")
  trt_lev <- levels(trt_fac)
  K <- length(trt_lev)
  if (K < 2L) stop("treatment column must have at least 2 levels.")

  cov_x <- stats::model.matrix(~ . - 1, data = data[, cov_vars, drop = FALSE])
  keep_balance <- which(apply(cov_x, 2, function(v) stats::var(v) > 0))
  if (length(keep_balance) == 0L) stop("All covariate balance columns have zero variance.")
  cov_x <- cov_x[, keep_balance, drop = FALSE]

  if (match_on == "gps") {
    loggps_names <- paste0("loggps_", trt_lev[-1L])
    missing_loggps <- setdiff(loggps_names, names(data))
    if (length(missing_loggps) > 0L) {
      stop(
        "Missing log-GPS columns: ",
        paste(missing_loggps, collapse = ", "),
        ". Set fit_gps=TRUE or pass data that already contains these columns."
      )
    }
    Z_all <- as.matrix(data[, loggps_names, drop = FALSE])
    if (anyNA(Z_all) || any(!is.finite(Z_all))) {
      stop("Found missing or non-finite loggps_* values. Set fit_gps=TRUE or check the supplied GPS columns.")
    }
    Z_all <- scale(Z_all)
    Z_all[is.na(Z_all)] <- 0
  } else {
    Z_all <- stats::model.matrix(~ . - 1, data = data[, cov_vars, drop = FALSE])
    keep <- which(apply(Z_all, 2, function(v) stats::var(v) > 0))
    if (length(keep) == 0L) stop("All covariate matching columns have zero variance.")
    Z_all <- Z_all[, keep, drop = FALSE]

    if (isTRUE(standardize)) {
      Z_all <- scale(Z_all)
      Z_all[is.na(Z_all)] <- 0
    }

    if (cov_distance == "mahalanobis") {
      S <- stats::cov(Z_all)
      scale0 <- mean(diag(S))
      if (!is.finite(scale0) || scale0 <= 0) scale0 <- 1
      diag(S) <- diag(S) + ridge * scale0
      L <- tryCatch(chol(solve(S)), error = function(e) NULL)
      if (is.null(L)) {
        eig <- eigen(S, symmetric = TRUE)
        L <- eig$vectors %*% diag(1 / sqrt(pmax(eig$values, ridge))) %*% t(eig$vectors)
      }
      Z_all <- Z_all %*% L
    }
  }

  idx_by_t <- lapply(trt_lev, function(tt) which(trt_fac == tt))
  names(idx_by_t) <- trt_lev

  matches <- vector("list", K)
  names(matches) <- trt_lev
  n <- nrow(data)

  for (k in seq_len(K)) {
    donors <- idx_by_t[[k]]
    if (length(donors) < match_ratio) {
      stop(sprintf("Not enough donors in arm '%s' for match_ratio=%d.", trt_lev[k], match_ratio))
    }
    not_t <- setdiff(seq_len(n), donors)
    mat <- matrix(NA_integer_, nrow = n, ncol = match_ratio)

    if (length(not_t) > 0L) {
      nn <- RANN::nn2(
        data = Z_all[donors, , drop = FALSE],
        query = Z_all[not_t, , drop = FALSE],
        k = match_ratio,
        treetype = "kd",
        searchtype = "priority"
      )$nn.idx
      if (match_ratio == 1L) nn <- matrix(nn, nrow = length(not_t), ncol = 1L)
      mat[not_t, ] <- matrix(donors[as.vector(nn)], nrow = length(not_t), ncol = match_ratio)
    }
    mat[donors, 1L] <- donors
    matches[[k]] <- mat
  }

  .smd <- function(x1, x2) {
    x1 <- as.numeric(x1)
    x2 <- as.numeric(x2)
    v1 <- stats::var(x1)
    v2 <- stats::var(x2)
    denom <- sqrt((v1 + v2) / 2)
    if (!is.finite(denom) || denom <= 0) return(0)
    abs(mean(x1) - mean(x2)) / denom
  }

  contrast_pairs <- utils::combn(trt_lev, 2, simplify = FALSE)
  rows <- list()
  row_id <- 0L

  for (pair in contrast_pairs) {
    a <- pair[1L]
    b <- pair[2L]
    ia <- idx_by_t[[a]]
    ib <- idx_by_t[[b]]
    contrast_name <- paste0(a, " vs ", b)

    matched_b_for_a <- as.vector(matches[[b]][ia, , drop = FALSE])
    matched_a_for_b <- as.vector(matches[[a]][ib, , drop = FALSE])
    matched_b_for_a <- matched_b_for_a[!is.na(matched_b_for_a)]
    matched_a_for_b <- matched_a_for_b[!is.na(matched_a_for_b)]

    ia_rep <- rep(ia, each = match_ratio)
    ib_rep <- rep(ib, each = match_ratio)
    ia_rep <- ia_rep[seq_along(matched_b_for_a)]
    ib_rep <- ib_rep[seq_along(matched_a_for_b)]

    for (j in seq_len(ncol(cov_x))) {
      cov_name <- colnames(cov_x)[j]
      row_id <- row_id + 1L
      rows[[row_id]] <- data.frame(
        contrast = contrast_name,
        covariate = cov_name,
        sample = "Before matching",
        abs_smd = .smd(cov_x[ia, j], cov_x[ib, j]),
        stringsAsFactors = FALSE
      )

      row_id <- row_id + 1L
      rows[[row_id]] <- data.frame(
        contrast = contrast_name,
        covariate = cov_name,
        sample = "After matching",
        abs_smd = .smd(
          c(cov_x[ia_rep, j], cov_x[matched_a_for_b, j]),
          c(cov_x[matched_b_for_a, j], cov_x[ib_rep, j])
        ),
        stringsAsFactors = FALSE
      )
    }
  }

  balance <- do.call(rbind, rows)
  balance$sample <- factor(balance$sample, levels = c("Before matching", "After matching"))

  if (style == "love") {
    plot_data <- stats::aggregate(
      abs_smd ~ covariate + sample,
      data = balance,
      FUN = max
    )
    plot_data$sample <- factor(
      ifelse(plot_data$sample == "Before matching", "All", "Matched"),
      levels = c("All", "Matched")
    )

    cov_order <- plot_data[plot_data$sample == "All", c("covariate", "abs_smd")]
    cov_order <- cov_order[order(cov_order$abs_smd, decreasing = FALSE), ]
    plot_data$covariate <- factor(plot_data$covariate, levels = cov_order$covariate)

    p <- ggplot2::ggplot(
      plot_data,
      ggplot2::aes(x = .data$abs_smd, y = .data$covariate,
                   shape = .data$sample, fill = .data$sample)
    ) +
      ggplot2::geom_vline(xintercept = 0, color = "black", linewidth = 0.35) +
      ggplot2::geom_vline(xintercept = threshold / 2, linetype = "dashed",
                          color = "black", linewidth = 0.35) +
      ggplot2::geom_vline(xintercept = threshold, color = "black", linewidth = 0.35) +
      ggplot2::geom_point(size = 2.8, color = "black", stroke = 0.8) +
      ggplot2::scale_shape_manual(values = c("All" = 21, "Matched" = 16)) +
      ggplot2::scale_fill_manual(values = c("All" = "white", "Matched" = "black")) +
      ggplot2::labs(
        x = "Absolute Standardized\nMean Difference",
        y = NULL,
        shape = NULL,
        fill = NULL
      ) +
      ggplot2::theme_bw() +
      ggplot2::theme(
        legend.position = c(0.84, 0.14),
        legend.background = ggplot2::element_rect(fill = "white", color = "black"),
        panel.grid.major.y = ggplot2::element_line(linetype = "dotted", color = "grey72"),
        panel.grid.major.x = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank()
      )
  } else {
    plot_data <- balance
    cov_order <- stats::aggregate(abs_smd ~ covariate, data = plot_data, FUN = max)
    cov_order <- cov_order[order(cov_order$abs_smd, decreasing = FALSE), ]
    plot_data$covariate <- factor(plot_data$covariate, levels = cov_order$covariate)

    p <- ggplot2::ggplot(
      plot_data,
      ggplot2::aes(x = .data$abs_smd, y = .data$covariate,
                   color = .data$sample, shape = .data$sample)
    ) +
      ggplot2::geom_vline(xintercept = threshold, linetype = "dashed", color = "grey55") +
      ggplot2::geom_point(
        position = ggplot2::position_dodge(width = 0.45),
        size = 2.2
      ) +
      ggplot2::facet_wrap(~contrast) +
      ggplot2::labs(
        x = "Absolute standardized mean difference",
        y = NULL,
        color = NULL,
        shape = NULL
      ) +
      ggplot2::theme_bw() +
      ggplot2::theme(
        legend.position = "bottom",
        panel.grid.major.y = ggplot2::element_line(color = "grey90"),
        panel.grid.minor = ggplot2::element_blank()
      )
  }

  if (isTRUE(return_data)) {
    return(list(plot = p, balance = balance, plot_data = plot_data))
  }

  p
}
