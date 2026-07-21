#' GPS- or Covariate-matching ATE estimator (numeric-column interface)
#'
#' For cross-fitted inference, use dr_gpsm(); gps_matching()'s bootstrap is intended for standalone (non-cross-fitted) usage.
#'
#' @param data        Data frame
#' @param treatment   Column index of the treatment variable (factor)
#' @param outcome     Column index of the outcome variable
#' @param pred        pred K × n matrix of m_k(X) predictions aligned with rows of data
#' @param contrast    C(K, 2) × K contrast matrix from build_contrast()
#' @param nboot       Number of bootstrap replications (only used if do_boot=TRUE; can be NULL otherwise)
#' @param match_on    'gps' (default) or 'covariates'
#' @param covariate   Numeric vector of covariate column indices (required if match_on='covariates')
#' @param cov_distance 'euclidean' or 'mahalanobis' (only for match_on='covariates')
#' @param standardize Logical; standardize covariate features before matching
#' @param ridge       Small ridge for covariance in Mahalanobis whitening
#' @param match_ratio Integer ≥ 1; number of matches per unit per target group
#' @param return_tau  Logical; if TRUE, also return tau (and tau_centered)
#' @param do_boot     Logical; if TRUE, compute bootstrap CI within this fold; if FALSE, skip (ci_* = NA)
#' @param boot_weight Bootstrap reweighting scheme: 'multinom' or 'exp' (only used when \code{do_boot=TRUE}).
#'
#' @return list(estimate, ci_lower, ci_upper, tau?, tau_centered?)
#' @export
#'
#' @examples
#' set.seed(1)
#' n <- 75
#' dat <- data.frame(
#'   trt = factor(rep(c("A", "B", "C"), each = 25)),
#'   x1 = rnorm(n),
#'   x2 = runif(n)
#' )
#' dat$y <- 1 + 0.4 * (dat$trt == "B") + 0.8 * (dat$trt == "C") +
#'   dat$x1 + rnorm(n, sd = 0.5)
#'
#' gps_dat <- gps_pre_process(
#'   dat,
#'   treatment = 1,
#'   treatment_ref = "A",
#'   covariate = 2:3,
#'   gps_model = "logit"
#' )
#' contrast <- build_contrast(levels(gps_dat$trt), ref = "A")
#' arm_means <- tapply(gps_dat$y, gps_dat$trt, mean)[levels(gps_dat$trt)]
#' pred <- matrix(rep(as.numeric(arm_means), each = nrow(gps_dat)),
#'   nrow = length(arm_means),
#'   byrow = TRUE
#' )
#'
#' gps_matching(
#'   data = gps_dat,
#'   treatment = 1,
#'   outcome = 4,
#'   pred = pred,
#'   contrast = contrast,
#'   do_boot = FALSE
#' )
gps_matching <- function(data,
                         treatment,
                         outcome,
                         pred,
                         contrast,
                         nboot = NULL,
                         match_on = c("gps","covariates"),
                         covariate = NULL,
                         cov_distance = c("mahalanobis","euclidean"),
                         standardize = FALSE,
                         ridge = 1e-8,
                         match_ratio = 1L,
                         return_tau = FALSE,
                         do_boot = FALSE,
                         boot_weight = c("multinom","exp")) {

  if (!requireNamespace("RANN", quietly = TRUE)) {
    stop("Package 'RANN' is required for matching (RANN::nn2).")
  }

  match_on     <- match.arg(match_on)
  cov_distance <- match.arg(cov_distance)
  boot_weight  <- match.arg(boot_weight)

  stopifnot(is.numeric(match_ratio), length(match_ratio) == 1L, match_ratio >= 1)
  match_ratio <- as.integer(match_ratio)

  if (isTRUE(do_boot)) {
    if (is.null(nboot) || length(nboot) != 1L || !is.finite(nboot) || nboot < 1) {
      stop("When do_boot=TRUE, please provide a valid nboot (integer >= 1).")
    }
    nboot <- as.integer(nboot)
  } else {
    nboot <- 0L
  }

  n <- nrow(data)
  stopifnot(is.factor(data[[treatment]]))
  tlev <- levels(data[[treatment]])
  K <- length(tlev)

  ## pred sanity: K x n
  pred <- as.matrix(pred)
  if (nrow(pred) != K) stop("`pred` must have nrow = number of treatment levels K.")
  if (ncol(pred) != n) stop("`pred` must have ncol equal to nrow(data) within this fold.")

  ## outcome
  y <- data[[outcome]]

  ## --- Build matching feature matrix Z_all (n x d) ---
  if (match_on == "gps") {
    loggps_cols <- integer(K - 1L)
    for (k in 2:K) {
      nm <- paste0("loggps_", tlev[k])
      idx <- match(nm, names(data))
      if (is.na(idx)) {
        stop(sprintf("Can't find %s (ensure gps_pre_process ran, and fold overwrote loggps_*).", nm))
      }
      loggps_cols[k - 1L] <- idx
    }

    Z_all <- as.matrix(data[, loggps_cols, drop = FALSE])
    Z_all <- scale(Z_all)
    Z_all[is.na(Z_all)] <- 0
    if (any(!is.finite(Z_all))) {
      stop("Found non-finite loggps_* in this fold. Did you overwrite fold-specific loggps_*?")
    }

  } else {
    if (is.null(covariate) || length(covariate) < 1L) {
      stop("When match_on='covariates', please provide `covariate` column indices.")
    }

    X <- stats::model.matrix(~ . - 1, data = data[, covariate, drop = FALSE])

    keep <- which(apply(X, 2, function(v) stats::var(v) > 0))
    if (length(keep) == 0L) stop("All covariate columns have zero variance in this fold.")
    X <- X[, keep, drop = FALSE]

    if (isTRUE(standardize)) {
      X <- scale(X)
      X[is.na(X)] <- 0
    }

    if (cov_distance == "mahalanobis") {
      S <- stats::cov(X)
      scale0 <- mean(diag(S))
      if (!is.finite(scale0) || scale0 <= 0) scale0 <- 1
      diag(S) <- diag(S) + ridge * scale0
      L <- tryCatch(chol(solve(S)), error = function(e) NULL)
      if (is.null(L)) {
        eig <- eigen(S, symmetric = TRUE)
        L <- eig$vectors %*% diag(1 / sqrt(pmax(eig$values, ridge))) %*% t(eig$vectors)
      }
      X <- X %*% L
    }

    Z_all <- as.matrix(X)
  }

  ## --- Indices by treatment ---
  T_fac <- data[[treatment]]
  idx_by_t <- lapply(tlev, function(tt) which(T_fac == tt))
  names(idx_by_t) <- tlev

  ## --- Build matches[[t]] for NOT-treated units only (treated units conceptually have J_t(i)={i}) ---
  matches <- vector("list", K)
  names(matches) <- tlev

  for (k in seq_len(K)) {
    donors  <- idx_by_t[[k]]
    treated <- donors

    if (length(donors) < match_ratio) {
      stop(sprintf("Not enough donors in arm '%s' for match_ratio=%d (only %d).",
                   tlev[k], match_ratio, length(donors)))
    }

    M <- min(match_ratio, length(donors))  # effective M in this fold/arm

    # We only need nearest neighbors for units with T != t
    not_t <- setdiff(seq_len(n), treated)

    # Initialize matches matrix (n x M) for consistency/debug; treated rows not used for omega counting
    mat <- matrix(NA_integer_, nrow = n, ncol = M)

    if (length(not_t) > 0) {
      Z_don <- Z_all[donors, , drop = FALSE]
      Z_q   <- Z_all[not_t, , drop = FALSE]

      nn <- RANN::nn2(data = Z_don, query = Z_q, k = M,
                      treetype = "kd", searchtype = "priority")$nn.idx
      if (M == 1L) nn <- matrix(nn, nrow = length(not_t), ncol = 1L)

      # Map donor-row indices back to global indices
      mat[not_t, ] <- matrix(donors[as.vector(nn)], nrow = length(not_t), ncol = M)
    }

    # Enforce algorithmic convention J_t(i) = {i} for treated i
    # (store i in first column; remaining columns left NA to avoid accidental multi-count)
    mat[treated, 1L] <- treated

    matches[[k]] <- mat
  }

  ## --- Compute psi_all (K x n) by the note formula (Step 4) ---
  psi_all <- matrix(NA_real_, nrow = K, ncol = n, dimnames = list(tlev, NULL))

  for (k in seq_len(K)) {
    treated <- idx_by_t[[k]]
    not_t   <- setdiff(seq_len(n), treated)
    M <- ncol(matches[[k]])  # effective M used in this arm

    # Count how often each i is used as a donor among NOT-treated units' matching sets
    match_count_not <- integer(n)
    if (length(not_t) > 0) {
      donor_mat <- matches[[k]][not_t, , drop = FALSE]
      donor_vec <- as.vector(donor_mat)
      donor_vec <- donor_vec[!is.na(donor_vec)]
      if (length(donor_vec) > 0) {
        match_count_not <- tabulate(donor_vec, nbins = n)
      }
    }

    # omega_i = (1/M) * sum_{ell in I_b} 1{i in J_t(ell)}
    # For treated i: J_t(i) = {i} contributes +1, so omega gets +1/M in addition to being used by not_t units.
    omega <- match_count_not / M

    coef_i <- 1 + omega
    a_t <- as.numeric(pred[k, ])

    # psi_{t,i} = a_{t,i} + 1(T_i=t)*(1+omega_i)*(Y_i - a_{t,i})
    psi_all[k, ] <- a_t
    if (length(treated) > 0) {
      psi_all[k, treated] <- a_t[treated] + coef_i[treated] * (y[treated] - a_t[treated])
    }
  }

  ## --- Contrast estimates and per-unit contributions ---
  P <- nrow(contrast)
  tau_contrib <- matrix(NA_real_, nrow = P, ncol = n, dimnames = list(rownames(contrast), NULL))
  tau_hat <- numeric(P); names(tau_hat) <- rownames(contrast)

  for (p in seq_len(P)) {
    neg <- which(contrast[p, ] == -1)
    pos <- which(contrast[p, ] ==  1)
    if (length(neg) != 1L || length(pos) != 1L) stop("Each contrast row must have exactly one -1 and one +1.")
    tau_contrib[p, ] <- psi_all[pos, ] - psi_all[neg, ]
    tau_hat[p] <- mean(tau_contrib[p, ])
  }

  ## --- Optional bootstrap CI ---
  if (isTRUE(do_boot)) {
    delta_star <- matrix(NA_real_, nrow = nboot, ncol = P)
    for (b in seq_len(nboot)) {
      W <- if (boot_weight == "exp") {
        stats::rexp(n, rate = 1)
      } else {
        as.numeric(stats::rmultinom(1, size = n, prob = rep(1/n, n)))
      }
      Wc <- W - mean(W)
      delta_star[b, ] <- drop(crossprod(Wc, t(tau_contrib))) / n
    }
    q_lo <- apply(delta_star, 2, stats::quantile, probs = 0.025, names = FALSE)
    q_hi <- apply(delta_star, 2, stats::quantile, probs = 0.975, names = FALSE)
    ci_lower <- tau_hat - q_hi
    ci_upper <- tau_hat - q_lo
  } else {
    ci_lower <- rep(NA_real_, P)
    ci_upper <- rep(NA_real_, P)
  }

  out <- list(
    estimate = stats::setNames(as.numeric(tau_hat), rownames(contrast)),
    ci_lower = stats::setNames(as.numeric(ci_lower), rownames(contrast)),
    ci_upper = stats::setNames(as.numeric(ci_upper), rownames(contrast))
  )

  if (isTRUE(return_tau)) {
    out$tau <- tau_contrib
    out$psi <- psi_all
  }

  out
}
