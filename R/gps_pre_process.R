#' Generalized Propensity Score (GPS) Pre-processing (with optional GPS fitting)
#'
#' In the cross-fitted BC-GPSM workflow, this function should primarily:
#'  1) enforce treatment factor + reference level;
#'  2) keep only complete cases (treatment + covariates);
#'  3) ALWAYS create gps_ and loggps_ placeholder columns (NA by default);
#'  4) optionally fit a full-sample GPS model (diagnostic only, e.g., for histogram).
#'
#' @param data A data.frame containing treatment/covariates (and possibly outcome).
#' @param treatment Single numeric column index for treatment.
#' @param treatment_ref Optional reference level. If \code{NULL}, the last
#'   observed treatment level is used. The selected reference is moved to the
#'   first factor level internally because the log-GPS ratios use it as the
#'   denominator.
#' @param covariate Numeric vector of covariate column indices.
#' @param gps_model One of "logit","gbm","gam","xgboost","ranger"
#'   (used only if fit_gps=TRUE). Flexible learners require optional packages;
#'   the XGBoost and ranger backends are experimental and under active testing.
#' @param gps_params Optional named list of model-specific parameters.
#' @param tune Logical; if TRUE, optional \pkg{caret} tuning is used for
#'   "logit" and "gbm" when fit_gps=TRUE.
#' @param tune_control Optional \code{caret::trainControl()} object for tuning.
#' @param tune_grids Optional tuning grids list.
#' @param seed Integer seed.
#' @param fit_gps Logical; if FALSE, do NOT fit GPS; return placeholders only.
#'
#' @return data.frame with original columns + gps_ and loggps_ (always present).
#' @export
#' @importFrom stats na.omit relevel reformulate as.formula
#' @importFrom splines bs
#' @importFrom nnet multinom
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
#' gps_dat <- gps_pre_process(
#'   data = dat,
#'   treatment = 1,
#'   covariate = 2:3,
#'   gps_model = "logit"
#' )
#'
#' levels(gps_dat$trt)
#' head(gps_dat[, grep("^(gps|loggps)_", names(gps_dat))])
gps_pre_process <- function(data,
                            treatment,
                            treatment_ref = NULL,
                            covariate,
                            gps_model = c("logit","gbm","gam","xgboost","ranger"),
                            gps_params = NULL,
                            tune = FALSE,
                            tune_control = NULL,
                            tune_grids = NULL,
                            seed = 12345,
                            fit_gps = TRUE) {

  ## -- helpers --
  `%||%` <- function(a, b) if (is.null(a)) b else a

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

  get_user <- function(name) {
    if (is.null(gps_params) || !is.list(gps_params)) return(NULL)
    gps_params[[name]]
  }

  .default_trainControl <- function() {
    .drgpsm_require_optional("caret", "tune=TRUE")
    caret::trainControl(
      method = "cv", number = 5,
      classProbs = TRUE,
      summaryFunction = caret::multiClassSummary,
      savePredictions = "final",
      allowParallel = TRUE
    )
  }

  .default_grids <- function(model_name) {
    switch(model_name,
           "logit" = expand.grid(decay = 10^seq(-4, -1, length.out = 5)),
           "gbm"   = expand.grid(
             n.trees = c(1500, 3000, 4500),
             interaction.depth = c(2, 3, 4),
             shrinkage = c(0.01, 0.005),
             n.minobsinnode = c(5, 10)
           ),
           NULL
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
        'Invalid tuneGrid for method "%s"; expected columns: %s',
        method, paste(expected, collapse = ", ")
      ))
    }
    invisible(TRUE)
  }

  clip_prob <- function(p, eps = 1e-6) pmax(pmin(p, 1 - eps), eps)

  ## -- input checks --
  if (!is.data.frame(data)) stop("`data` must be a data.frame.")
  if (!is.numeric(treatment) || length(treatment) != 1L)
    stop("`treatment` must be a single numeric column index.")
  if (!is.numeric(covariate) || length(covariate) < 1L)
    stop("`covariate` must be a numeric vector of column indices (length >= 1).")

  p <- ncol(data)
  if (any(c(treatment, covariate) < 1 | c(treatment, covariate) > p)) {
    stop("treatment/covariate indices out of range.")
  }

  trt_var  <- names(data)[treatment]
  cov_vars <- names(data)[covariate]

  ## -- complete-case filtering for treatment + covariates (consistent with your estimator) --
  cc_cols <- unique(c(trt_var, cov_vars))
  data <- data[stats::complete.cases(data[, cc_cols, drop = FALSE]), , drop = FALSE]
  if (nrow(data) == 0L) stop("No complete cases after filtering on treatment + covariates.")

  ## -- enforce treatment factor + reference level --
  trt_fac <- factor(data[[trt_var]])
  original_levels <- levels(trt_fac)
  if (is.null(treatment_ref)) {
    treatment_ref <- original_levels[length(original_levels)]
  }
  if (!treatment_ref %in% original_levels) {
    stop("`treatment_ref` must be a level in treatment column.")
  }
  trt_fac <- stats::relevel(trt_fac, ref = treatment_ref)
  data[[trt_var]] <- trt_fac

  K <- nlevels(trt_fac)
  if (K < 2L) stop("treatment column must have at least 2 levels.")

  gps_model <- match.arg(gps_model)
  set.seed(seed)

  ## -- ALWAYS create placeholders (NA) --
  levs <- levels(trt_fac)
  gps_cols <- matrix(NA_real_, nrow = nrow(data), ncol = K,
                     dimnames = list(NULL, paste0("gps_", levs)))
  loggps_cols <- matrix(NA_real_, nrow = nrow(data), ncol = K - 1L,
                        dimnames = list(NULL, paste0("loggps_", levs[-1])))

  ## -- optionally fit GPS (diagnostic only, not used for cross-fitted estimation) --
  if (isTRUE(fit_gps)) {

    ## NOTE: This full-sample GPS is NOT cross-fitted. In your dr_gpsm you overwrite
    ## loggps_* on each fold's df_test using fold-specific PS.
    gps_raw <- NULL

    if (gps_model == "logit") {

      form <- stats::reformulate(cov_vars, response = trt_var)

      if (isTRUE(tune)) {
        tc   <- tune_control %||% .default_trainControl()
        grid <- if (!is.null(tune_grids) && !is.null(tune_grids$logit)) tune_grids$logit else .default_grids("logit")

        ## caret method name is "multinom"
        .validate_tunegrid("multinom", grid)

        user_par <- get_user("logit"); user_par <- if (is.null(user_par)) list() else user_par

        fit <- caret::train(
          form, data = data, method = "multinom",
          metric = "logLoss", trControl = tc, tuneGrid = grid,
          trace = FALSE,
          MaxNWts = user_par$MaxNWts %||% 10000,
          maxit   = user_par$maxit   %||% 200
        )

        gps_raw <- stats::predict(fit, newdata = data, type = "prob")
        gps_raw <- as.matrix(gps_raw)[, levs, drop = FALSE]

      } else {

        logit_def <- list(trace = FALSE)
        logit_par <- .merge_defaults(logit_def, get_user("logit"))

        fit  <- do.call(nnet::multinom, c(list(formula = form, data = data), logit_par))
        gps_raw <- stats::fitted(fit)

        ## ensure K-column matrix
        gps_raw <- as.matrix(gps_raw)
        if (ncol(gps_raw) != K) {
          # in rare cases, fitted() might not return all columns; force by predict(type="probs")
          gps_raw <- stats::predict(fit, newdata = data, type = "probs")
          gps_raw <- as.matrix(gps_raw)
        }
        colnames(gps_raw) <- levs
      }

    } else if (gps_model == "gbm") {

      .drgpsm_require_optional("gbm", "gps_model='gbm'")

      form <- stats::reformulate(cov_vars, response = trt_var)

      if (isTRUE(tune)) {
        tc   <- tune_control %||% .default_trainControl()
        grid <- if (!is.null(tune_grids) && !is.null(tune_grids$gbm)) tune_grids$gbm else .default_grids("gbm")

        ## caret method name is "gbm"
        .validate_tunegrid("gbm", grid)

        fit <- caret::train(
          form, data = data, method = "gbm",
          distribution = "multinomial",
          metric = "logLoss", trControl = tc, tuneGrid = grid,
          verbose = FALSE
        )

        gps_raw <- stats::predict(fit, newdata = data, type = "prob")
        gps_raw <- as.matrix(gps_raw)[, levs, drop = FALSE]

      } else {

        gbm_def <- list(
          n.trees = 3000L, interaction.depth = 3L, shrinkage = 0.01,
          cv.folds = 5L, n.minobsinnode = 10L, n.cores = 1L, verbose = FALSE
        )
        gbm_par <- .merge_defaults(gbm_def, get_user("gbm"))
        filt <- .filter_to_formals(gbm::gbm, gbm_par)

        fit <- suppressWarnings(do.call(
          gbm::gbm,
          c(list(formula = form, data = data, distribution = "multinomial"), filt$keep)
        ))

        best_iter <- if (!is.null(gbm_par$cv.folds) && gbm_par$cv.folds > 1L) {
          gbm::gbm.perf(fit, method = "cv", plot.it = FALSE)
        } else {
          gbm_par$n.trees
        }

        pred_arr <- gbm::predict.gbm(fit, newdata = data, n.trees = best_iter, type = "response")

        ## gbm multinomial often returns array: n x K x 1 OR n x K
        if (is.array(pred_arr) && length(dim(pred_arr)) >= 2L) {
          # try to coerce to n x K
          if (length(dim(pred_arr)) == 3L) {
            gps_raw <- pred_arr[, , 1, drop = FALSE][, , 1]
          } else {
            gps_raw <- pred_arr
          }
        } else {
          gps_raw <- pred_arr
        }
        gps_raw <- as.matrix(gps_raw)
        colnames(gps_raw) <- levs
      }

    } else if (gps_model == "xgboost") {

      if (isTRUE(tune)) {
        stop("GPS tuning is not currently implemented for gps_model='xgboost'; pass settings through gps_params$xgboost instead.")
      }

      gps_raw <- .drgpsm_fit_xgboost_gps(
        train_df = data,
        eval_df = data,
        t_name = trt_var,
        x_names = cov_vars,
        tr_levels = levs,
        params = get_user("xgboost"),
        seed = seed
      )

    } else if (gps_model == "ranger") {

      if (isTRUE(tune)) {
        stop("GPS tuning is not currently implemented for gps_model='ranger'; pass settings through gps_params$ranger instead.")
      }

      gps_raw <- .drgpsm_fit_ranger_gps(
        train_df = data,
        eval_df = data,
        t_name = trt_var,
        x_names = cov_vars,
        tr_levels = levs,
        params = get_user("ranger"),
        seed = seed
      )

    } else if (gps_model == "gam") {

      .drgpsm_require_optional("VGAM", "gps_model='gam' diagnostic fitting")

      make_term <- function(v) {
        x <- data[[v]]
        if (!is.numeric(x)) return(v)
        u <- length(unique(stats::na.omit(x)))
        if (u <= 3L) v else sprintf("splines::bs(%s, df=%d)", v, min(5L, u - 1L))
      }
      rhs  <- vapply(cov_vars, make_term, character(1))
      form <- stats::as.formula(paste(trt_var, "~", paste(rhs, collapse = " + ")))

      fit <- suppressWarnings(VGAM::vglm(
        formula = form,
        family  = VGAM::multinomial(),
        data    = data
      ))

      gps_raw <- stats::predict(fit, type = "response")
      gps_raw <- as.matrix(gps_raw)
      colnames(gps_raw) <- levs
    }

    ## finalize gps/loggps (diagnostic)
    gps_raw <- clip_prob(gps_raw, eps = 1e-6)
    gps_cols[,] <- gps_raw
    loggps_cols[,] <- log(gps_raw[, -1, drop = FALSE] / gps_raw[, 1])
  }

  ## bind + return
  out <- cbind(data, gps_cols, loggps_cols)
  out
}
