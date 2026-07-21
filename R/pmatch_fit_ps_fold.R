#' Fit PS model on training fold and predict probabilities on eval fold
#'
#' @param train_df training data frame
#' @param eval_df  evaluation data frame
#' @param t_name   treatment column name (factor)
#' @param x_names  covariate column names
#' @param model "logit","gbm","gam","xgboost","ranger". Flexible learners
#'   require optional packages; XGBoost and ranger are experimental.
#' @param tune Logical; optional caret tuning for logit/gbm.
#' @param tune_control \code{caret::trainControl()} or NULL.
#' @param tune_grid named list with components logit/gbm (same style as your gps_tune_grids)
#' @param gps_params optional named list of model params (same style as your gps_params)
#' @param tr_levels full set of treatment levels to enforce (character)
#' @param seed integer seed
#' @return list(e_hat_eval = n_eval x K matrix, levels = tr_levels)
#' @export
#'
#' @importFrom stats density
#' @importFrom utils modifyList
#'
#' @examples
#' set.seed(1)
#' n <- 75
#' dat <- data.frame(
#'   trt = factor(rep(c("A", "B", "C"), each = 25)),
#'   x1 = rnorm(n),
#'   x2 = runif(n)
#' )
#'
#' train_id <- c(1:20, 26:45, 51:70)
#' ps_fit <- pmatch_fit_ps_fold(
#'   train_df = dat[train_id, ],
#'   eval_df = dat[-train_id, ],
#'   t_name = "trt",
#'   x_names = c("x1", "x2"),
#'   model = "logit",
#'   tr_levels = levels(dat$trt)
#' )
#'
#' head(ps_fit$e_hat_eval)
pmatch_fit_ps_fold <- function(train_df, eval_df,
                               t_name, x_names,
                               model = c("logit","gbm","gam","xgboost","ranger"),
                               tune = FALSE,
                               tune_control = NULL,
                               tune_grid = NULL,
                               gps_params = NULL,
                               tr_levels = NULL,
                               seed = 12345) {
  model <- match.arg(model)

  `%||%` <- function(a,b) if (is.null(a)) b else a

  get_user <- function(name) {
    if (is.null(gps_params) || !is.list(gps_params)) return(NULL)
    gps_params[[name]]
  }

  # enforce levels
  if (is.null(tr_levels)) {
    train_df[[t_name]] <- factor(train_df[[t_name]])
    tr_levels <- levels(train_df[[t_name]])
  }
  train_df[[t_name]] <- factor(train_df[[t_name]], levels = tr_levels)
  eval_df[[t_name]]  <- factor(eval_df[[t_name]],  levels = tr_levels)

  fml <- stats::as.formula(paste(t_name, "~", paste(x_names, collapse = "+")))
  set.seed(seed)

  .multi_logloss <- function(data, lev = NULL, model = NULL) {
    y <- data$obs
    p <- as.matrix(data[, lev, drop = FALSE])
    c(logLoss = ModelMetrics::mlogLoss(y, p))
  }

  .default_tc <- function() {
    .drgpsm_require_optional("caret", "tune=TRUE")
    .drgpsm_require_optional("ModelMetrics", "GPS log-loss tuning")
    caret::trainControl(
      method="cv", number=5,
      classProbs=TRUE,
      summaryFunction=.multi_logloss,
      savePredictions="final",
      allowParallel=TRUE
    )
  }

  .validate_tunegrid <- function(method, grid) {
    .drgpsm_require_optional("caret", "tune=TRUE")
    mi_all <- caret::getModelInfo(method)
    if (length(mi_all) == 0L || is.null(mi_all[[method]])) {
      stop(sprintf("caret method '%s' not found in getModelInfo()", method))
    }
    expected <- mi_all[[method]]$parameters$parameter
    if (!is.data.frame(grid) || !setequal(colnames(grid), expected)) {
      stop(sprintf(
        'Invalid tuneGrid for method "%s". Expected columns: %s',
        method, paste(expected, collapse = ", ")
      ))
    }
    invisible(TRUE)
  }

  ## --- logit ---
  if (model == "logit") {
    if (isTRUE(tune)) {
      tc <- tune_control %||% .default_tc()
      grid <- NULL
      if (!is.null(tune_grid) && !is.null(tune_grid$logit)) grid <- tune_grid$logit
      if (is.null(grid)) grid <- expand.grid(decay = 10^seq(-4,-1,length.out=5))

      user_par <- get_user("logit"); if (is.null(user_par)) user_par <- list()

      .validate_tunegrid("multinom", grid)

      fit <- caret::train(
        fml, data = train_df, method = "multinom",
        metric = "logLoss",
        trControl = tc,
        tuneGrid  = grid,
        trace = FALSE,
        MaxNWts = user_par$MaxNWts %||% 10000,
        maxit   = user_par$maxit %||% 200
      )
      e_hat <- stats::predict(fit, newdata = eval_df, type = "prob")
      e_hat <- as.matrix(e_hat)
    } else {
      requireNamespace("nnet", quietly = TRUE)
      logit_def <- list(trace = FALSE)
      logit_par <- get_user("logit"); if (is.null(logit_par)) logit_par <- list()
      args <- c(list(formula = fml, data = train_df), logit_def, logit_par)

      fit <- suppressWarnings(do.call(nnet::multinom, args))
      e_hat <- stats::predict(fit, newdata = eval_df, type = "probs")

      if (is.null(dim(e_hat))) {
        e_hat <- cbind(1 - e_hat, e_hat)
        colnames(e_hat) <- tr_levels
      } else {
        e_hat <- e_hat[, tr_levels, drop = FALSE]
      }
    }

    e_hat <- as.matrix(e_hat)
    if (is.null(colnames(e_hat))) colnames(e_hat) <- tr_levels
    e_hat <- e_hat[, tr_levels, drop = FALSE]
    return(list(e_hat_eval = e_hat, levels = tr_levels))
  }

  ## --- gbm ---
  if (model == "gbm") {
    if (isTRUE(tune)) {
      tc <- tune_control %||% .default_tc()
      grid <- NULL
      if (!is.null(tune_grid) && !is.null(tune_grid$gbm)) grid <- tune_grid$gbm
      if (is.null(grid)) {
        grid <- expand.grid(
          n.trees = c(1500, 3000, 4500),
          interaction.depth = c(2,3,4),
          shrinkage = c(0.01, 0.005),
          n.minobsinnode = c(5,10)
        )
      }

      .validate_tunegrid("gbm", grid)
      fit <- caret::train(
        fml, data = train_df, method = "gbm",
        distribution = "multinomial",
        metric = "logLoss",
        trControl = tc,
        tuneGrid  = grid,
        verbose   = FALSE
      )
      e_hat <- stats::predict(fit, newdata = eval_df, type = "prob")
      e_hat <- as.matrix(e_hat)[, tr_levels, drop = FALSE]
      return(list(e_hat_eval = e_hat, levels = tr_levels))
    }

    .drgpsm_require_optional("gbm", "model='gbm'")
    gbm_def <- list(
      n.trees = 5000L,
      interaction.depth = 2L,
      shrinkage = 0.005,
      bag.fraction = 0.5,
      cv.folds = 5L,
      n.minobsinnode = 20L,
      n.cores = 1L,
      verbose = FALSE
    )

    gbm_par <- get_user("gbm"); if (is.null(gbm_par)) gbm_par <- list()
    par <- modifyList(gbm_def, gbm_par)

    fit <- suppressWarnings(gbm::gbm(
      formula = fml, data = train_df,
      distribution = "multinomial",
      n.trees = par$n.trees,
      interaction.depth = par$interaction.depth,
      shrinkage = par$shrinkage,
      cv.folds = par$cv.folds,
      n.minobsinnode = par$n.minobsinnode,
      n.cores = par$n.cores,
      verbose = par$verbose,
      bag.fraction = par$bag.fraction
    ))

    best_iter <- if (!is.null(par$cv.folds) && par$cv.folds > 1L) {
      gbm::gbm.perf(fit, method = "cv",  plot.it = FALSE)
    } else {
      gbm::gbm.perf(fit, method = "OOB", plot.it = FALSE)
    }
    best_iter <- max(1L, as.integer(best_iter))

    pred <- gbm::predict.gbm(
      fit,
      newdata = eval_df,
      n.trees = best_iter,
      type = "response"
    )

    if (is.array(pred) && length(dim(pred)) == 3L) {
      pred <- pred[, , dim(pred)[3], drop = FALSE]
      pred <- pred[, , 1]
    }

    e_hat <- as.matrix(pred)

    gbm_classes <- fit$classes
    if (is.null(gbm_classes)) {
      gbm_classes <- colnames(e_hat)
    }

    if (is.null(gbm_classes) || length(gbm_classes) != length(tr_levels)) {
      stop("GBM PS: cannot determine class order (fit$classes/colnames missing or wrong length).")
    }

    colnames(e_hat) <- gbm_classes

    if (!setequal(colnames(e_hat), tr_levels)) {
      stop(sprintf(
        "GBM PS: class set mismatch.\n  gbm: %s\n  expected: %s",
        paste(colnames(e_hat), collapse = ", "),
        paste(tr_levels, collapse = ", ")
      ))
    }

    e_hat <- e_hat[, tr_levels, drop = FALSE]

    rs <- rowSums(e_hat)
    if (any(!is.finite(rs)) || any(rs <= 0)) stop("GBM PS: non-finite/zero row sums in predicted probs.")

    eps <- 1e-12
    e_hat <- pmax(e_hat, eps)
    e_hat <- e_hat / rowSums(e_hat)

    return(list(e_hat_eval = e_hat, levels = tr_levels))
  }

  ## --- xgboost ---
  if (model == "xgboost") {
    if (isTRUE(tune)) {
      stop("GPS tuning is not currently implemented for model='xgboost'; pass settings through gps_params$xgboost instead.")
    }
    e_hat <- .drgpsm_fit_xgboost_gps(
      train_df = train_df,
      eval_df = eval_df,
      t_name = t_name,
      x_names = x_names,
      tr_levels = tr_levels,
      params = get_user("xgboost"),
      seed = seed
    )
    return(list(e_hat_eval = e_hat, levels = tr_levels))
  }

  ## --- ranger ---
  if (model == "ranger") {
    if (isTRUE(tune)) {
      stop("GPS tuning is not currently implemented for model='ranger'; pass settings through gps_params$ranger instead.")
    }
    e_hat <- .drgpsm_fit_ranger_gps(
      train_df = train_df,
      eval_df = eval_df,
      t_name = t_name,
      x_names = x_names,
      tr_levels = tr_levels,
      params = get_user("ranger"),
      seed = seed
    )
    return(list(e_hat_eval = e_hat, levels = tr_levels))
  }

  ## --- gam ---
  if (model == "gam") {
    .drgpsm_require_optional("mgcv", "model='gam'")

    # enforce levels
    train_df[[t_name]] <- factor(train_df[[t_name]], levels = tr_levels)
    eval_df[[t_name]]  <- factor(eval_df[[t_name]],  levels = tr_levels)

    K <- length(tr_levels)
    if (K < 2L) stop("Need at least 2 treatment levels for GAM PS.")

    eps <- 1e-6
    ref <- tr_levels[1]

    # ----- standardize numeric X using TRAIN stats -----
    num_vars <- x_names[vapply(train_df[, x_names, drop = FALSE], is.numeric, logical(1))]
    scale_par <- lapply(num_vars, function(v) {
      m <- mean(train_df[[v]], na.rm = TRUE)
      s <- stats::sd(train_df[[v]], na.rm = TRUE)
      if (!is.finite(s) || s <= 0) s <- 1
      list(mean = m, sd = s)
    })
    names(scale_par) <- num_vars

    for (v in num_vars) {
      m <- scale_par[[v]]$mean
      s <- scale_par[[v]]$sd
      train_df[[v]] <- (train_df[[v]] - m) / s
      eval_df[[v]]  <- (eval_df[[v]]  - m) / s
      train_df[[v]][!is.finite(train_df[[v]])] <- 0
      eval_df[[v]][!is.finite(eval_df[[v]])]   <- 0
    }

    # ----- smooth specification -----
    gam_user <- get_user("gam"); if (is.null(gam_user)) gam_user <- list()
    df_max <- gam_user$df_max %||% 6L
    min_unique_smooth <- gam_user$min_unique_smooth %||% 5L
    method_use <- gam_user$method %||% "REML"
    maxit_use  <- gam_user$maxit %||% 200L

    make_term <- function(v, df_for_k) {
      x <- df_for_k[[v]]
      if (!is.numeric(x)) return(v)
      u <- length(unique(stats::na.omit(x)))
      if (u < min_unique_smooth) return(v)
      k_use <- min(as.integer(df_max), as.integer(u - 1L))
      k_use <- max(3L, k_use)
      sprintf("s(%s, k=%d)", v, k_use)
    }

    # ----- fit K-1 binomial GAMs for log(e_k/e_ref) -----
    eta <- matrix(0, nrow = nrow(eval_df), ncol = K - 1L)
    colnames(eta) <- tr_levels[-1]

    for (kk in 2:K) {
      levk <- tr_levels[kk]

      df_k <- train_df[train_df[[t_name]] %in% c(ref, levk), , drop = FALSE]
      if (nrow(df_k) < 10L) stop(sprintf("Too few samples for GAM PS: ref=%s vs %s", ref, levk))

      # binary response: 1 for levk, 0 for ref
      df_k[[".ybin"]] <- as.integer(df_k[[t_name]] == levk)

      rhs <- vapply(x_names, function(v) make_term(v, df_k), character(1))
      form <- stats::as.formula(paste(".ybin ~", paste(rhs, collapse = " + ")))

      fit_k <- mgcv::gam(
        formula = form,
        family  = stats::binomial(),
        data    = df_k,
        method  = method_use,
        control = mgcv::gam.control(maxit = as.integer(maxit_use))
      )

      # log-odds for levk vs ref
      eta[, kk - 1L] <- as.numeric(stats::predict(fit_k, newdata = eval_df, type = "link"))
    }

    # ----- softmax to probabilities -----
    # e_ref = 1 / (1 + sum exp(eta_k))
    # e_k   = exp(eta_k) * e_ref
    exp_eta <- exp(pmin(pmax(eta, -30), 30))
    denom <- 1 + rowSums(exp_eta)
    e_ref <- 1 / denom

    P <- matrix(NA_real_, nrow = nrow(eval_df), ncol = K)
    colnames(P) <- tr_levels
    P[, 1] <- e_ref
    for (kk in 2:K) {
      P[, kk] <- exp_eta[, kk - 1L] * e_ref
    }

    # clamp + renorm
    P[!is.finite(P)] <- eps
    P <- pmin(pmax(P, eps), 1 - eps)
    P <- P / rowSums(P)

    return(list(e_hat_eval = P, levels = tr_levels))
  }

  stop("Unknown PS model.")
}
