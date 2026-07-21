#' Doubly-robust GPS matching estimator
#'
#' Implements a doubly-robust generalized propensity score (GPS) matching estimator
#' with optional cross-fitting and bootstrap inference. GPS can be estimated via
#' multinomial logistic regression, GBM, GAM, XGBoost, or ranger; outcome models
#' can be fitted via LM, RF, GBM, GAM, XGBoost, or ranger. Flexible learners and
#' tuning integrations require optional packages. The XGBoost and ranger
#' backends are experimental and remain under active testing.
#'
#' @details
#' The treatment, covariate, and outcome arguments must identify distinct,
#' named columns. Missing outcome values are not currently modeled and produce
#' an error with guidance to remove or impute them before fitting. For
#' `outcome_model = "lm"`, binary outcomes should be encoded numerically as
#' 0/1. Each treatment arm must contain at least `folds` observations and
#' at least `match_ratio` potential donors. Invalid settings are reported
#' before nuisance-model fitting with an actionable error message.
#'
#' @param data Data frame with treatment, covariates, and outcome variables
#' @param treatment Column index of the treatment variable
#' @param treatment_ref Reference treatment level. If \code{NULL}, the last
#'   observed treatment level is used. The selected reference is moved to the
#'   first factor level internally for GPS log-ratio construction.
#' @param covariate Numeric vector of covariate column indices
#' @param outcome Column index of the outcome variable
#' @param gps_model Choice of GPS model: 'logit', 'gbm', 'gam', 'xgboost', or
#'   'ranger'. All choices except 'logit' require optional packages.
#' @param outcome_model Choice of outcome model: 'none','lm','rf','gbm','gam',
#'   'xgboost', or 'ranger'. Flexible learners require optional packages.
#' @param folds Number of folds for cross-fitting (default = 2)
#' @param nboot Number of bootstrap replications (default = 500)
#' @param hist Logical; TRUE to save a diagnostic full-sample GPS histogram
#' @param hist_path File path (pdf/png/jpg) to save histogram if hist = TRUE
#' @param gps_params Optional named list of hyperparameters per GPS model
#' @param outcome_params Optional named list of hyperparameters per outcome model
#' @param match_on 'gps' (default) or 'covariates'
#' @param cov_distance For covariate matching: 'euclidean' or 'mahalanobis' (default 'mahalanobis')
#' @param standardize If TRUE, center/scale covariate features before matching (default FALSE)
#' @param match_ratio Integer >= 1; number of matches per unit per target group (default 1)
#' @param gps_tune Logical; if TRUE, enable optional caret tuning for GPS
#'   ('logit','gbm').
#' @param gps_tune_control Optional \code{caret::trainControl()} for GPS tuning; default is 5-fold CV
#' @param gps_tune_grids Optional named list of GPS tuning grids (e.g., \code{list(gbm = ..., logit = ...)})
#' @param gps_seed Integer seed for GPS tuning and fitting (default 12345)
#' @param outcome_tune Logical; if TRUE, enable optional caret tuning for outcome
#'   models ('rf','gbm').
#' @param outcome_tune_control Optional \code{caret::trainControl()} for outcome tuning
#' @param outcome_tune_grids Optional named list of outcome tuning grids (e.g., \code{list(rf = ..., gbm = ...)})
#' @param outcome_seed Integer seed for outcome tuning and fitting (default 12345)
#' @param fold_seed Integer seed for fold assignment (default 12345).
#' @param boot_weight Bootstrap reweighting scheme: 'multinom' or 'exp'.
#' @param two_step_calibration Logical; if TRUE, apply two-step residual calibration using ridge on GPS index V.
#' @param calib_shrinkage Numeric in \code{[0,1]}; shrinkage multiplier applied to the calibration step size gamma.
#'
#' @return A list with components:
#' \describe{
#'   \item{estimate}{Named numeric vector of doubly-robust ATE estimates for each contrast.}
#'   \item{ci_lower}{Named numeric vector of lower 95% confidence bounds.}
#'   \item{ci_upper}{Named numeric vector of upper 95% confidence bounds.}
#' }
#'
#' @export
#' @importFrom stats na.omit
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
#' fit <- dr_gpsm(
#'   data = dat,
#'   treatment = 1,
#'   treatment_ref = "A",
#'   covariate = 2:3,
#'   outcome = 4,
#'   gps_model = "logit",
#'   outcome_model = "lm",
#'   folds = 2,
#'   nboot = 10
#' )
#'
#' fit$estimate
dr_gpsm <- function(data,
                    treatment,
                    treatment_ref = NULL,
                    covariate,
                    outcome,
                    gps_model = c("logit","gbm","gam","xgboost","ranger"),
                    outcome_model = c("none","lm","rf","gbm","gam","xgboost","ranger"),
                    folds = 2,
                    fold_seed = 12345,
                    nboot = 500,
                    boot_weight = c("multinom","exp"),
                    hist = FALSE,
                    hist_path = NULL,
                    gps_params = NULL,
                    outcome_params = NULL,
                    match_on = c("gps","covariates"),
                    cov_distance = c("mahalanobis","euclidean"),
                    standardize = FALSE,
                    match_ratio = 1L,
                    gps_tune = FALSE,
                    gps_tune_control = NULL,
                    gps_tune_grids = NULL,
                    gps_seed = 12345,
                    outcome_tune = FALSE,
                    outcome_tune_control = NULL,
                    outcome_tune_grids = NULL,
                    outcome_seed = 12345,
                    two_step_calibration = FALSE,
                    calib_shrinkage = 1) {

  ## -- helpers --
  .merge_defaults <- function(defaults, user) {
    if (is.null(user)) return(defaults)
    stopifnot(is.list(user))
    for (nm in names(user)) defaults[[nm]] <- user[[nm]]
    defaults
  }
  .filter_to_formals <- function(fun, args) {
    keep <- intersect(names(args), names(formals(fun)))
    list(keep = args[keep], dropped = setdiff(names(args), keep))
  }
  get_outcome_user <- function(name) {
    if (is.null(outcome_params) || !is.list(outcome_params)) return(NULL)
    outcome_params[[name]]
  }
  `%||%` <- function(a,b) if (is.null(a)) b else a

  # separate caret defaults for classification vs regression
  .default_tc_cls <- function(summary = c("multi","two")) {
    .drgpsm_require_optional("caret", "model tuning")
    summary <- match.arg(summary)
    caret::trainControl(
      method = "cv", number = 5,
      classProbs = TRUE,
      summaryFunction = switch(summary,
                               multi = caret::multiClassSummary,
                               two   = caret::twoClassSummary
      ),
      savePredictions = "final",
      allowParallel   = FALSE
    )
  }

  .default_tc_reg <- function() {
    .drgpsm_require_optional("caret", "model tuning")
    caret::trainControl(
      method = "cv", number = 5,
      classProbs = FALSE,
      summaryFunction = caret::defaultSummary, # RMSE/RSquared/MAE
      savePredictions = "final",
      allowParallel   = TRUE
    )
  }

  .validate_tunegrid <- function(method, grid) {
    .drgpsm_require_optional("caret", "model tuning")
    mi_all <- caret::getModelInfo(method)
    if (length(mi_all) == 0L || is.null(mi_all[[method]])) {
      stop(sprintf("caret method '%s' not found in getModelInfo()", method))
    }
    expected <- mi_all[[method]]$parameters$parameter
    if (!is.data.frame(grid) || !setequal(colnames(grid), expected)) {
      stop(sprintf('Invalid tuneGrid for method "%s"; please use getModelInfo("%s")[[1]]$parameters to build the proper tuning grid',
                   method, method))
    }
    invisible(TRUE)
  }

  # coercion: enforce binary outcome as factor(one/zero), with "one" as event (first level)
  .coerce_binary_outcome <- function(df, outcome_var) {
    yy <- df[[outcome_var]]
    yy_chr <- as.character(yy)

    low <- tolower(yy_chr)
    yy_chr <- ifelse(low %in% c("1","yes","true","case","event","one"), "one",
                     ifelse(low %in% c("0","no","false","control","non-event","nonevent","zero"), "zero", yy_chr))

    if (!all(stats::na.omit(yy_chr) %in% c("one","zero"))) {
      yy_num <- suppressWarnings(as.numeric(yy_chr))
      if (all(stats::na.omit(yy_num) %in% c(0,1))) {
        yy_chr <- ifelse(yy_num == 1, "one", "zero")
      } else {
        stop("Binary outcome detected, but cannot coerce outcome values to {0,1} or {one,zero}.")
      }
    }

    df[[outcome_var]] <- factor(yy_chr, levels = c("one","zero"))
    df
  }

  .coerce_binary_outcome01 <- function(df, outcome_var) {
    yy <- df[[outcome_var]]
    yy_chr <- as.character(yy)
    low <- tolower(yy_chr)
    yy_chr <- ifelse(low %in% c("1","yes","true","case","event","one"), "1",
                     ifelse(low %in% c("0","no","false","control","non-event","nonevent","zero"), "0", yy_chr))
    yy_num <- suppressWarnings(as.numeric(yy_chr))
    if (!all(stats::na.omit(yy_num) %in% c(0,1))) {
      stop("Binary outcome detected, but cannot coerce to numeric {0,1} for gbm.")
    }
    df[[outcome_var]] <- yy_num
    df
  }

  .safe_prob_one <- function(prob_mat) {
    prob_mat <- as.data.frame(prob_mat)

    if ("one" %in% colnames(prob_mat)) {
      return(as.numeric(prob_mat[["one"]]))
    }

    if (ncol(prob_mat) == 2L) {
      # fallback: prefer column 2 (often event prob); but "one" should exist if we coerced correctly
      return(as.numeric(prob_mat[[2L]]))
    }

    stop("predict(..., type='prob') returned no column named 'one' and not a 2-class matrix.")
  }


  .y_to_01 <- function(yvec, is_binary) {
    if (!is_binary) return(as.numeric(yvec))
    if (is.factor(yvec)) {
      # you coerced binary outcome to levels c("one","zero") in your tuning helpers
      return(as.numeric(yvec == levels(yvec)[1L]))  # first level as event
    }
    # numeric 0/1 already
    as.numeric(yvec)
  }

  .is_whole_number <- function(x) {
    is.numeric(x) && length(x) == 1L && !is.na(x) && is.finite(x) &&
      x == floor(x)
  }

  .validate_column_index <- function(index, argument, multiple = FALSE) {
    valid_length <- if (multiple) length(index) >= 1L else length(index) == 1L
    if (!is.numeric(index) || !valid_length || anyNA(index) ||
        any(!is.finite(index)) || any(index != floor(index))) {
      expected <- if (multiple) "one or more integer column indices" else "one integer column index"
      stop(sprintf("`%s` must contain %s.", argument, expected), call. = FALSE)
    }
    if (any(index < 1L | index > ncol(data))) {
      stop(
        sprintf(
          "`%s` contains a column index outside 1:%d. Check the column positions in `data`.",
          argument, ncol(data)
        ),
        call. = FALSE
      )
    }
    invisible(TRUE)
  }

  ## -- args and preprocessing --
  if (requireNamespace("foreach", quietly = TRUE)) foreach::registerDoSEQ()

  if (!is.data.frame(data)) {
    stop("`data` must be a data.frame.", call. = FALSE)
  }
  if (nrow(data) == 0L || ncol(data) == 0L) {
    stop("`data` must contain at least one row and one column.", call. = FALSE)
  }
  if (is.null(names(data)) || any(!nzchar(names(data))) || anyDuplicated(names(data))) {
    stop("Every column in `data` must have a unique, non-empty name.", call. = FALSE)
  }
  .validate_column_index(treatment, "treatment")
  .validate_column_index(covariate, "covariate", multiple = TRUE)
  .validate_column_index(outcome, "outcome")
  if (anyDuplicated(covariate)) {
    stop("`covariate` must not contain duplicate column indices.", call. = FALSE)
  }
  if (treatment %in% covariate || outcome %in% covariate || treatment == outcome) {
    stop(
      "`treatment`, `covariate`, and `outcome` must refer to distinct columns.",
      call. = FALSE
    )
  }

  gps_model     <- match.arg(gps_model)
  outcome_model <- match.arg(outcome_model)
  match_on      <- match.arg(match_on)
  cov_distance  <- match.arg(cov_distance)

  if (!.is_whole_number(folds) || folds < 2L) {
    stop("`folds` must be a single integer greater than or equal to 2.", call. = FALSE)
  }
  folds <- as.integer(folds)
  if (!.is_whole_number(nboot) || nboot < 1L) {
    stop("`nboot` must be a single integer greater than or equal to 1.", call. = FALSE)
  }
  nboot <- as.integer(nboot)
  if (!.is_whole_number(match_ratio) || match_ratio < 1L) {
    stop("`match_ratio` must be a single integer greater than or equal to 1.", call. = FALSE)
  }
  match_ratio <- as.integer(match_ratio)

  boot_weight <- match.arg(boot_weight)
  if (!is.logical(two_step_calibration) || length(two_step_calibration) != 1L ||
      is.na(two_step_calibration)) {
    stop("`two_step_calibration` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.numeric(calib_shrinkage) || length(calib_shrinkage) != 1L ||
      is.na(calib_shrinkage) || !is.finite(calib_shrinkage) ||
      calib_shrinkage < 0 || calib_shrinkage > 1) {
    stop("`calib_shrinkage` must be one number between 0 and 1.", call. = FALSE)
  }

  outcome_var_input <- names(data)[outcome]
  missing_outcomes <- sum(is.na(data[[outcome_var_input]]))
  if (missing_outcomes > 0L) {
    suffix <- if (missing_outcomes == 1L) "value" else "values"
    stop(
      sprintf(
        "The `outcome` column '%s' contains %d missing %s. `dr_gpsm()` does not currently model missing outcomes; remove or impute them before fitting.",
        outcome_var_input, missing_outcomes, suffix
      ),
      call. = FALSE
    )
  }
  if (outcome_model == "lm" && !is.numeric(data[[outcome_var_input]])) {
    stop(
      "`outcome_model = 'lm'` requires a numeric outcome. Encode a binary outcome as numeric 0/1, or select a model that supports factor outcomes.",
      call. = FALSE
    )
  }

  if (isTRUE(two_step_calibration) && outcome_model == "none") {
    stop(
      "`two_step_calibration = TRUE` requires an outcome model; choose `outcome_model = 'lm'`, 'rf', 'gbm', 'gam', 'xgboost', or 'ranger'.",
      call. = FALSE
    )
  }
  if (isTRUE(two_step_calibration)) {
    .drgpsm_require_optional("glmnet", "two_step_calibration=TRUE")
  }

  gps_tune_control_use <- gps_tune_control
  if (isTRUE(gps_tune) && is.null(gps_tune_control_use)) {
    gps_tune_control_use <- .default_tc_cls("multi")
  }

  dat_processed <- gps_pre_process(
    data = data,
    treatment = treatment,
    treatment_ref = treatment_ref,
    covariate = covariate,
    gps_model = gps_model,
    gps_params = gps_params,
    tune = gps_tune,
    tune_control = gps_tune_control_use,
    tune_grids = gps_tune_grids,
    seed = gps_seed,
    fit_gps = FALSE
  )

  if (isTRUE(hist)) {
    if (is.null(hist_path)) stop("When hist=TRUE, please provide hist_path.")
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
      stop("Package 'ggplot2' is required when hist=TRUE.")
    }
    dat_hist <- gps_pre_process(
      data = data,
      treatment = treatment,
      treatment_ref = treatment_ref,
      covariate = covariate,
      gps_model = gps_model,
      gps_params = gps_params,
      tune = gps_tune,
      tune_control = gps_tune_control_use,
      tune_grids = gps_tune_grids,
      seed = gps_seed,
      fit_gps = TRUE
    )
    p <- gps_histogram(dat_hist)
    ggplot2::ggsave(filename = hist_path, plot = p, width = 8, height = 8)
  }

  trt_var  <- names(data)[treatment]
  trt_lev  <- levels(dat_processed[[trt_var]])
  if (is.null(treatment_ref)) treatment_ref <- trt_lev[1L]
  contrast <- build_contrast(trt_lev, ref = treatment_ref)

  # ensure reference level aligns with trt_lev[1]
  if (trt_lev[1] != treatment_ref) {
    stop(sprintf(
      "Internal: treatment_ref='%s' is not the first level after preprocessing (first='%s').",
      treatment_ref, trt_lev[1]
    ))
  }

  n <- nrow(dat_processed)
  arm_counts <- table(dat_processed[[trt_var]])
  if (any(arm_counts < folds)) {
    bad_arm <- names(arm_counts)[which(arm_counts < folds)[1L]]
    bad_n <- unname(arm_counts[[bad_arm]])
    suffix <- if (bad_n == 1L) "observation" else "observations"
    stop(
      sprintf(
        "Cannot create reliable cross-fitting folds: treatment arm '%s' has %d %s, but `folds = %d`. Each treatment arm needs at least `folds` observations. Reduce `folds`, combine or remove the rare arm, or use more data.",
        bad_arm, bad_n, suffix, folds
      ),
      call. = FALSE
    )
  }
  if (any(arm_counts < match_ratio)) {
    bad_arm <- names(arm_counts)[which(arm_counts < match_ratio)[1L]]
    bad_n <- unname(arm_counts[[bad_arm]])
    stop(
      sprintf(
        "Treatment arm '%s' has %d observations, fewer than `match_ratio = %d`. Reduce `match_ratio` or use more observations in that arm.",
        bad_arm, bad_n, match_ratio
      ),
      call. = FALSE
    )
  }

  make_folds <- function(fac, K) {
    stopifnot(is.factor(fac))
    inds  <- split(seq_along(fac), fac)
    folds <- vector("list", K)
    for (ids in inds) {
      ids <- sample(ids)
      sp  <- split(ids, rep(1:K, length.out = length(ids)))
      for (f in seq_len(K)) folds[[f]] <- c(folds[[f]], sp[[f]])
    }
    lapply(folds, sort)
  }
  trt_fac  <- dat_processed[[trt_var]]
  set.seed(fold_seed)
  folds_id <- make_folds(trt_fac, folds)

  ## -- map original variable names to column indices in dat_processed --
  cov_vars    <- names(data)[covariate]
  outcome_var <- names(data)[outcome]

  treat_col_dp   <- match(trt_var, names(dat_processed))
  outcome_col_dp <- match(outcome_var, names(dat_processed))
  cov_cols_dp    <- match(cov_vars, names(dat_processed))

  if (is.na(treat_col_dp))   stop("Internal error: treatment column not found in dat_processed.")
  if (is.na(outcome_col_dp)) stop("Internal error: outcome column not found in dat_processed.")
  if (any(is.na(cov_cols_dp))) stop("Internal error: some covariate columns not found in dat_processed.")

  y <- dat_processed[[outcome_var]]
  is_binary <- {
    if (is.factor(y)) nlevels(y) == 2L else {
      u <- unique(stats::na.omit(y)); length(u) == 2L && all(sort(u) %in% c(0,1))
    }
  }

  ## -- containers for full-sample OOF nuisances --
  K <- length(trt_lev)

  # OOF GPS index used for matching: n x (K-1)
  loggps_oof <- matrix(NA_real_, nrow = n, ncol = K - 1L)
  colnames(loggps_oof) <- paste0("loggps_", trt_lev[-1L])

  # OOF outcome regression predictions: K x n
  pred_oof <- matrix(NA_real_, nrow = K, ncol = n,
                     dimnames = list(trt_lev, NULL))

  ## -- cross-fitting loop --
  for (f in seq_along(folds_id)) {
    test  <- folds_id[[f]]
    train <- setdiff(seq_len(n), test)

    ## ---- fold-specific OOF GPS / outcome nuisance estimation on eval fold ----
    df_train <- dat_processed[train, , drop = FALSE]
    df_test  <- dat_processed[test,  , drop = FALSE]

    tab_tr <- table(factor(df_train[[trt_var]], levels = trt_lev))
    if (any(tab_tr == 0L)) stop(sprintf("Fold %d: some treatment has 0 training samples.", f))

    ps_out <- pmatch_fit_ps_fold(
      train_df = df_train,
      eval_df  = df_test,
      t_name   = trt_var,
      x_names  = cov_vars,
      model    = gps_model,
      tune     = gps_tune,
      tune_control = gps_tune_control_use,
      tune_grid    = gps_tune_grids,
      gps_params   = gps_params,
      tr_levels    = trt_lev,
      seed         = gps_seed + f
    )
    e_hat_eval <- ps_out$e_hat_eval

    # normalize + build log-ratio relative to reference level (trt_lev[1])
    eps_ps <- 1e-6
    e_hat_eval <- pmax(e_hat_eval, eps_ps)
    e_hat_eval <- e_hat_eval / rowSums(e_hat_eval)

    loggps_raw <- log(e_hat_eval[, -1, drop = FALSE] / e_hat_eval[, 1])

    # for matching: cap
    log_cap <- 10
    loggps_match <- pmin(pmax(loggps_raw, -log_cap), log_cap)

    # for two-step: use raw V (NOT capped)
    V_hat_eval <- if (isTRUE(two_step_calibration) && outcome_model != "none") as.matrix(loggps_raw) else NULL

    # Always define pred_test (K x n_test), even if outcome_model == "none"
    pred_test <- matrix(0, nrow = length(trt_lev), ncol = nrow(df_test))

    if (outcome_model == "none") {
      y_train_num_all <- .y_to_01(df_train[[outcome_var]], is_binary)
      for (k in seq_along(trt_lev)) {
        idx_k <- which(df_train[[trt_var]] == trt_lev[k])
        mk <- mean(y_train_num_all[idx_k])
        pred_test[k, ] <- mk
      }
    } else {

      set.seed(outcome_seed + f)

      fit_fun <-
        switch(outcome_model,
               lm = {
                 function(df) stats::lm(stats::reformulate(cov_vars, response = outcome_var), data = df)
               },

               rf = {
                 if (isTRUE(outcome_tune)) {

                   tc <- if (!is_binary) outcome_tune_control %||% .default_tc_reg()
                   else                  outcome_tune_control %||% .default_tc_cls("two")

                   .drgpsm_require_optional("caret", "outcome_tune=TRUE")
                   .drgpsm_require_optional("randomForest", "outcome_model='rf'")
                   if (is_binary) .drgpsm_require_optional("pROC", "binary outcome tuning")

                   grid <- NULL
                   if (!is.null(outcome_tune_grids) && !is.null(outcome_tune_grids$rf)) {
                     grid <- outcome_tune_grids$rf
                     .validate_tunegrid("rf", grid)
                   } else {
                     grid <- expand.grid(mtry = pmax(1, floor(sqrt(length(cov_vars))) + c(-1, 0, 1)))
                     .validate_tunegrid("rf", grid)
                   }

                   function(df) {
                     df2 <- df
                     if (is_binary) df2 <- .coerce_binary_outcome(df2, outcome_var)
                     form <- stats::reformulate(cov_vars, response = outcome_var)

                     caret::train(
                       form, data = df2,
                       method    = "rf",
                       trControl = tc,
                       tuneGrid  = grid,
                       metric    = if (is_binary) "ROC" else "RMSE"
                     )
                   }

                 } else {

                   rf_def <- rf_def <- list(
                     ntree = 2000L,
                     mtry = NULL,
                     nodesize = 5L,
                     sampsize = NULL
                   )
                   rf_par <- .merge_defaults(rf_def, get_outcome_user("rf"))

                   function(df) {
                     .drgpsm_require_optional("randomForest", "outcome_model='rf'")

                     df2 <- df
                     if (is_binary) df2 <- .coerce_binary_outcome(df2, outcome_var)

                     form <- stats::reformulate(cov_vars, response = outcome_var)

                     p <- length(cov_vars)
                     n <- nrow(df2)

                     if (is.null(rf_par$mtry)) {
                       rf_par$mtry <- if (is_binary) max(1L, floor(sqrt(p))) else max(1L, floor(p / 3))
                     }

                     args <- c(list(formula = form, data = df2), rf_par)

                     args <- args[!vapply(args, is.null, logical(1))]

                     do.call(randomForest::randomForest, args)
                   }
                 }
               },

               gbm = {
                 if (isTRUE(outcome_tune)) {

                   tc <- if (!is_binary) outcome_tune_control %||% .default_tc_reg()
                   else                  outcome_tune_control %||% .default_tc_cls("two")

                   .drgpsm_require_optional("caret", "outcome_tune=TRUE")
                   .drgpsm_require_optional("gbm", "outcome_model='gbm'")
                   if (is_binary) .drgpsm_require_optional("pROC", "binary outcome tuning")

                   grid <- NULL
                   if (!is.null(outcome_tune_grids) && !is.null(outcome_tune_grids$gbm)) {
                     grid <- outcome_tune_grids$gbm
                     .validate_tunegrid("gbm", grid)
                   } else {
                     grid <- expand.grid(
                       n.trees = c(1500, 3000),
                       interaction.depth = c(2, 3, 4),
                       shrinkage = c(0.01, 0.005),
                       n.minobsinnode = c(5, 10)
                     )
                     .validate_tunegrid("gbm", grid)
                   }

                   function(df) {
                     df2 <- df
                     if (is_binary) df2 <- .coerce_binary_outcome(df2, outcome_var)

                     form <- stats::reformulate(cov_vars, response = outcome_var)

                     caret::train(
                       form, data = df2,
                       method       = "gbm",
                       distribution = if (is_binary) "bernoulli" else "gaussian",
                       metric       = if (is_binary) "ROC" else "RMSE",
                       trControl    = tc,
                       tuneGrid     = grid,
                       verbose      = FALSE
                     )
                   }

                 } else {

                   gbm_def <- list(
                     n.trees = 3000L,
                     interaction.depth = 3L,
                     shrinkage = 0.01,
                     cv.folds = 5L,
                     n.minobsinnode = 10L,
                     verbose = FALSE
                   )
                   gbm_par <- .merge_defaults(gbm_def, get_outcome_user("gbm"))

                   function(df) {
                     .drgpsm_require_optional("gbm", "outcome_model='gbm'")

                     df2 <- df
                     if (is_binary) df2 <- .coerce_binary_outcome01(df2, outcome_var)

                     form  <- stats::reformulate(cov_vars, response = outcome_var)
                     distn <- if (is_binary) "bernoulli" else "gaussian"

                     filt <- .filter_to_formals(gbm::gbm, gbm_par)

                     do.call(
                       gbm::gbm,
                       c(list(formula = form, data = df2, distribution = distn), filt$keep)
                     )
                   }
                 }
               },

               xgboost = {
                 if (isTRUE(outcome_tune)) {
                   stop("Outcome tuning is not currently implemented for outcome_model='xgboost'; pass settings through outcome_params$xgboost instead.")
                 }
                 xgb_par <- .merge_defaults(
                   .drgpsm_xgboost_outcome_defaults(is_binary, outcome_seed + f),
                   get_outcome_user("xgboost")
                 )

                 function(df) {
                   y_train <- .y_to_01(df[[outcome_var]], is_binary)
                   .drgpsm_fit_xgboost_outcome(
                     df = df,
                     x_names = cov_vars,
                     y = y_train,
                     is_binary = is_binary,
                     params = xgb_par,
                     seed = outcome_seed + f
                   )
                 }
               },

               ranger = {
                 if (isTRUE(outcome_tune)) {
                   stop("Outcome tuning is not currently implemented for outcome_model='ranger'; pass settings through outcome_params$ranger instead.")
                 }
                 ranger_par <- .merge_defaults(
                   .drgpsm_ranger_outcome_defaults(length(cov_vars), is_binary, outcome_seed + f),
                   get_outcome_user("ranger")
                 )

                 function(df) {
                   df2 <- df
                   if (is_binary) df2 <- .coerce_binary_outcome(df2, outcome_var)
                   .drgpsm_fit_ranger_outcome(
                     df = df2,
                     x_names = cov_vars,
                     outcome_var = outcome_var,
                     is_binary = is_binary,
                     params = ranger_par,
                     seed = outcome_seed + f
                   )
                 }
               },

               gam = {
                 .drgpsm_require_optional("mgcv", "outcome_model='gam'")
                 gam_def <- list(df_max = 6L, family = NULL, method = "REML", min_unique_smooth = 5L)
                 gam_par <- .merge_defaults(gam_def, get_outcome_user("gam"))

                 fam_use <- if (!is.null(gam_par$family)) gam_par$family else
                   if (is_binary) stats::binomial() else stats::gaussian()

                 function(df) {
                   rhs_terms <- vapply(cov_vars, function(v) {
                     x <- df[[v]]
                     if (!is.numeric(x)) return(v)
                     u <- length(unique(stats::na.omit(x)))
                     if (u < gam_par$min_unique_smooth) return(v)
                     k_use <- min(as.integer(gam_par$df_max), as.integer(u - 1L))
                     k_use <- max(3L, k_use)
                     sprintf("s(%s, k=%d)", v, k_use)
                   }, character(1))

                   form <- stats::as.formula(paste(outcome_var, "~", paste(rhs_terms, collapse = " + ")))
                   mgcv::gam(formula = form, family = fam_use, data = df, method = gam_par$method)
                 }
               }
        )

      # 1) Fit within-treatment outcome models; predict on df_test only
      for (k in seq_along(trt_lev)) {

        df_k <- df_train[df_train[[trt_var]] == trt_lev[k], , drop = FALSE]
        if (nrow(df_k) == 0L) {
          stop(sprintf("Fold %d: no training observations for treatment '%s'. Increase folds or check sample size per arm.",
                       f, trt_lev[k]))
        }

        fit_k <- fit_fun(df_k)

        # ---- predict on df_test (existing) ----
        if (outcome_model %in% c("rf","gbm") && isTRUE(outcome_tune)) {

          if (is_binary) {
            prob_mat <- stats::predict(fit_k, newdata = df_test, type = "prob")
            p <- .safe_prob_one(prob_mat)
          } else {
            p <- as.numeric(stats::predict(fit_k, newdata = df_test, type = "raw"))
          }
          pred_test[k, ] <- as.numeric(p)

        } else if (outcome_model == "gbm") {

          ntrees <- if (!is.null(fit_k$cv.folds) && fit_k$cv.folds > 0L) {
            gbm::gbm.perf(fit_k, method = "cv", plot.it = FALSE)
          } else fit_k$n.trees

          pred_test[k, ] <- as.numeric(gbm::predict.gbm(fit_k, newdata = df_test,
                                                        n.trees = ntrees, type = "response"))

        } else if (outcome_model == "xgboost") {

          pred_test[k, ] <- .drgpsm_predict_xgboost_outcome(
            fit_k,
            newdata = df_test,
            x_names = cov_vars
          )

        } else if (outcome_model == "lm") {

          pred_test[k, ] <- as.numeric(stats::predict(fit_k, newdata = df_test, type = "response"))

        } else if (outcome_model == "rf") {

          if (is_binary) {
            prob_mat <- stats::predict(fit_k, newdata = df_test, type = "prob")
            pred_test[k, ] <- as.numeric(.safe_prob_one(prob_mat))

          } else {
            pred_test[k, ] <- as.numeric(stats::predict(fit_k, newdata = df_test, type = "response"))
          }

        } else if (outcome_model == "ranger") {

          pred_test[k, ] <- .drgpsm_predict_ranger_outcome(fit_k, newdata = df_test)

        } else if (outcome_model == "gam") {

          pred_test[k, ] <- as.numeric(stats::predict(fit_k, newdata = df_test, type = "response"))
        }
      }

      # 2) Two-step residual calibration (optional, note-style):
      #    For each arm t, fit ridge on eval-fold treated residuals R = Y - m_t(X),
      #    using V_hat (GPS index) as features, then update a_{t,i} = m_t(X_i) + s * r_t(V_i)
      #    for ALL i in the eval fold (s = calib_shrinkage).
      #
      #    Update: default lambda choice = lambda.1se (more stable),
      #            calib_shrinkage default handled by function signature (now 1),
      #            allow fixed lambda via outcome_params$two_step$lambda (numeric).
      if (isTRUE(two_step_calibration)) {

        y_test_num <- .y_to_01(df_test[[outcome_var]], is_binary)

        if (is.null(V_hat_eval)) stop("Internal error: V_hat_eval not prepared.")

        # ---- read optional two-step options from outcome_params (no interface change) ----
        two_step_user <- NULL
        if (!is.null(outcome_params) && is.list(outcome_params)) {
          # accept a few possible keys for robustness
          if (!is.null(outcome_params$two_step)) two_step_user <- outcome_params$two_step
          if (is.null(two_step_user) && !is.null(outcome_params$calibration)) two_step_user <- outcome_params$calibration
          if (is.null(two_step_user) && !is.null(outcome_params$two_step_calibration)) two_step_user <- outcome_params$two_step_calibration
        }
        if (!is.null(two_step_user)) stopifnot(is.list(two_step_user))

        # optional: fixed lambda numeric, or rule "lambda.1se"/"lambda.min"
        lambda_fixed <- if (!is.null(two_step_user$lambda)) two_step_user$lambda else NULL
        lambda_rule  <- if (!is.null(two_step_user$lambda_rule)) two_step_user$lambda_rule else "lambda.1se"
        lambda_rule  <- match.arg(lambda_rule, choices = c("lambda.1se","lambda.min"))

        # optional CV controls
        cv_nfolds <- if (!is.null(two_step_user$nfolds)) as.integer(two_step_user$nfolds) else 5L
        if (!is.finite(cv_nfolds) || cv_nfolds < 2L) cv_nfolds <- 5L

        for (k in seq_along(trt_lev)) {

          idx_eval_t <- which(df_test[[trt_var]] == trt_lev[k])
          # too few treated in eval fold -> skip calibration for this arm
          if (length(idx_eval_t) < 10L) next

          if (is_binary) pred_test[k, ] <- pmin(1, pmax(0, pred_test[k, ]))
          R_eval   <- y_test_num[idx_eval_t] - pred_test[k, idx_eval_t]
          V_eval_t <- V_hat_eval[idx_eval_t, , drop = FALSE]

          # drop columns with ~0 variance to avoid numerical issues (esp. small folds)
          vv <- apply(V_eval_t, 2, stats::var)
          keepj <- which(is.finite(vv) & vv > 0)
          if (length(keepj) == 0L) next
          V_eval_t2 <- V_eval_t[, keepj, drop = FALSE]
          V_all2    <- V_hat_eval[, keepj, drop = FALSE]

          # Fit ridge calibration model
          if (!is.null(lambda_fixed)) {
            # fixed lambda: no CV, deterministic
            lam <- as.numeric(lambda_fixed)
            if (length(lam) != 1L || !is.finite(lam) || lam <= 0) {
              stop("two_step: outcome_params$two_step$lambda must be a single positive number.")
            }

            fit <- glmnet::glmnet(
              x = V_eval_t2, y = R_eval,
              alpha = 0, intercept = TRUE, standardize = TRUE,
              lambda = lam
            )
            r_hat_all <- as.numeric(stats::predict(fit, newx = V_all2, s = lam))

          } else {
            # CV ridge: choose lambda.1se by default (more stable)
            nf <- min(cv_nfolds, length(R_eval))
            nf <- max(2L, nf)
            cv <- glmnet::cv.glmnet(
              x = V_eval_t2, y = R_eval,
              alpha = 0, intercept = TRUE, standardize = TRUE,
              nfolds = nf)

            s_use <- switch(lambda_rule,
                            "lambda.1se" = cv$lambda.1se,
                            "lambda.min" = cv$lambda.min)

            r_hat_all <- as.numeric(stats::predict(cv, newx = V_all2, s = s_use))
          }

          # Update predictions in eval fold: a_t <- m_t + shrinkage * r_hat(V)
          pred_test[k, ] <- pred_test[k, ] + calib_shrinkage * r_hat_all

          # keep probabilities valid for binary outcomes
          if (is_binary) pred_test[k, ] <- pmin(1, pmax(0, pred_test[k, ]))
        }
      }
    }


    stopifnot(nrow(pred_test) == length(trt_lev))
    stopifnot(ncol(pred_test) == length(test))

    # store OOF GPS index for this eval fold
    loggps_oof[test, ] <- loggps_match

    # store OOF outcome predictions for this eval fold
    pred_oof[, test] <- pred_test
  }

  ## -- check full-sample OOF nuisances are complete --
  if (anyNA(loggps_oof)) {
    stop("Internal error: loggps_oof has NA entries; some test indices were not filled.")
  }
  if (anyNA(pred_oof)) {
    stop("Internal error: pred_oof has NA entries; some test indices were not filled.")
  }

  ## -- overwrite full-sample loggps_* with OOF GPS index --
  dat_match <- dat_processed
  for (k in 2:length(trt_lev)) {
    nm <- paste0("loggps_", trt_lev[k])
    dat_match[[nm]] <- as.numeric(loggps_oof[, k - 1L])
  }

  ## -- one full-sample matching run using OOF nuisances --
  g <- gps_matching(
    data = dat_match,
    treatment = treat_col_dp,
    outcome = outcome_col_dp,
    pred = pred_oof,
    contrast = contrast,
    match_on = match_on,
    covariate = cov_cols_dp,
    cov_distance = cov_distance,
    standardize = standardize,
    match_ratio = match_ratio,
    return_tau = TRUE,
    do_boot = FALSE
  )

  stopifnot(
    is.matrix(g$tau),
    nrow(g$tau) == nrow(contrast),
    ncol(g$tau) == n
  )

  tau_all <- g$tau
  est <- rowMeans(tau_all)

  # Note-style linear-form bootstrap on per-unit contributions
  P <- nrow(contrast)
  delta_star <- matrix(NA_real_, nrow = nboot, ncol = P)

  # tau_all is P x n; transpose once for fast crossprod
  tau_t <- t(tau_all)  # n x P

  for (b in seq_len(nboot)) {
    W <- if (boot_weight == "exp") {
      stats::rexp(n, rate = 1)
    } else {
      as.numeric(stats::rmultinom(1, size = n, prob = rep(1/n, n)))
    }
    Wc <- W - mean(W)

    # delta* = (1/n) sum_i Wc_i * tau_i  (vector length P)
    delta_star[b, ] <- drop(crossprod(Wc, tau_t)) / n
  }

  q_lo <- apply(delta_star, 2, stats::quantile, probs = 0.025, names = FALSE)
  q_hi <- apply(delta_star, 2, stats::quantile, probs = 0.975, names = FALSE)

  ci_lower <- est - q_hi
  ci_upper <- est - q_lo

  list(
    estimate = stats::setNames(as.numeric(est), rownames(contrast)),
    ci_lower = stats::setNames(as.numeric(ci_lower), rownames(contrast)),
    ci_upper = stats::setNames(as.numeric(ci_upper), rownames(contrast))
  )
}
