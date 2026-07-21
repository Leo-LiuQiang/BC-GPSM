.drgpsm_require_optional <- function(package, feature) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop(
      sprintf(
        "Package '%s' is required for %s. Install this optional dependency with install.packages('%s').",
        package, feature, package
      ),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

.drgpsm_null_default <- function(a, b) if (is.null(a)) b else a
`%||%` <- .drgpsm_null_default

.drgpsm_model_matrix <- function(data, x_names, columns = NULL) {
  form <- stats::reformulate(x_names)
  x <- stats::model.matrix(form, data = data)
  x <- x[, colnames(x) != "(Intercept)", drop = FALSE]

  if (!is.null(columns)) {
    missing_cols <- setdiff(columns, colnames(x))
    if (length(missing_cols) > 0L) {
      zeros <- matrix(0, nrow = nrow(x), ncol = length(missing_cols),
                      dimnames = list(NULL, missing_cols))
      x <- cbind(x, zeros)
    }
    x <- x[, columns, drop = FALSE]
  }

  if (ncol(x) == 0L) {
    stop("No usable covariate columns after constructing the model matrix.")
  }
  storage.mode(x) <- "double"
  x[!is.finite(x)] <- 0
  x
}

.drgpsm_prob_matrix <- function(prob, levels, label, eps = 1e-12) {
  prob <- as.matrix(prob)
  if (ncol(prob) != length(levels)) {
    stop(sprintf(
      "%s: expected %d probability columns, got %d.",
      label, length(levels), ncol(prob)
    ))
  }
  if (is.null(colnames(prob))) colnames(prob) <- levels
  if (!setequal(colnames(prob), levels)) {
    stop(sprintf(
      "%s: probability class mismatch.\n  model: %s\n  expected: %s",
      label,
      paste(colnames(prob), collapse = ", "),
      paste(levels, collapse = ", ")
    ))
  }
  prob <- prob[, levels, drop = FALSE]
  prob[!is.finite(prob)] <- eps
  prob <- pmax(prob, eps)
  prob <- prob / rowSums(prob)
  prob
}

.drgpsm_prob_one <- function(prob_mat) {
  prob_mat <- as.data.frame(prob_mat)
  if ("one" %in% colnames(prob_mat)) return(as.numeric(prob_mat[["one"]]))
  if (ncol(prob_mat) == 2L) return(as.numeric(prob_mat[[2L]]))
  stop("Predicted probability matrix has no column named 'one' and is not two-class.")
}

.drgpsm_xgboost_gps_defaults <- function(seed) {
  list(
    nrounds = 500L,
    max_depth = 3L,
    eta = 0.03,
    min_child_weight = 2,
    subsample = 0.85,
    colsample_bytree = 0.85,
    lambda = 1,
    alpha = 0,
    gamma = 0,
    tree_method = "hist",
    nthread = 1L,
    seed = seed,
    verbosity = 0,
    verbose = 0
  )
}

.drgpsm_xgboost_outcome_defaults <- function(is_binary, seed) {
  list(
    nrounds = if (is_binary) 450L else 550L,
    max_depth = 3L,
    eta = 0.03,
    min_child_weight = if (is_binary) 2 else 3,
    subsample = 0.85,
    colsample_bytree = 0.85,
    lambda = 1,
    alpha = 0,
    gamma = 0,
    tree_method = "hist",
    nthread = 1L,
    seed = seed,
    verbosity = 0,
    verbose = 0
  )
}

.drgpsm_xgboost_train <- function(x, y, params, objective, eval_metric, num_class = NULL) {
  .drgpsm_require_optional("xgboost", "model='xgboost'")

  nrounds <- as.integer(params$nrounds)
  verbose <- as.integer(params$verbose %||% 0)
  xgb_params <- params[setdiff(names(params), c("nrounds", "verbose"))]
  xgb_params$objective <- objective
  xgb_params$eval_metric <- eval_metric
  if (!is.null(num_class)) xgb_params$num_class <- num_class

  dtrain <- xgboost::xgb.DMatrix(data = x, label = y)
  xgboost::xgb.train(
    params = xgb_params,
    data = dtrain,
    nrounds = nrounds,
    verbose = verbose
  )
}

.drgpsm_fit_xgboost_gps <- function(train_df, eval_df, t_name, x_names,
                                    tr_levels, params = NULL, seed = 12345) {
  .drgpsm_require_optional("xgboost", "gps_model='xgboost'")

  train_x <- .drgpsm_model_matrix(train_df, x_names)
  eval_x <- .drgpsm_model_matrix(eval_df, x_names, columns = colnames(train_x))
  y <- as.integer(factor(train_df[[t_name]], levels = tr_levels)) - 1L

  defaults <- .drgpsm_xgboost_gps_defaults(seed)
  params <- utils::modifyList(defaults, params %||% list())
  fit <- .drgpsm_xgboost_train(
    x = train_x,
    y = y,
    params = params,
    objective = "multi:softprob",
    eval_metric = "mlogloss",
    num_class = length(tr_levels)
  )

  pred <- stats::predict(fit, newdata = xgboost::xgb.DMatrix(data = eval_x))
  pred <- matrix(pred, ncol = length(tr_levels), byrow = TRUE)
  colnames(pred) <- tr_levels
  .drgpsm_prob_matrix(pred, tr_levels, "XGBoost GPS")
}

.drgpsm_ranger_gps_defaults <- function(p, seed) {
  list(
    num.trees = 1000L,
    mtry = max(1L, floor(sqrt(p))),
    min.node.size = 10L,
    sample.fraction = 0.85,
    replace = FALSE,
    respect.unordered.factors = "order",
    oob.error = FALSE,
    num.threads = 1L,
    seed = seed,
    verbose = FALSE
  )
}

.drgpsm_fit_ranger_gps <- function(train_df, eval_df, t_name, x_names,
                                   tr_levels, params = NULL, seed = 12345) {
  .drgpsm_require_optional("ranger", "gps_model='ranger'")

  defaults <- .drgpsm_ranger_gps_defaults(length(x_names), seed)
  params <- utils::modifyList(defaults, params %||% list())
  params$probability <- NULL

  form <- stats::reformulate(x_names, response = t_name)
  args <- c(
    list(formula = form, data = train_df, probability = TRUE),
    params
  )
  args <- args[!vapply(args, is.null, logical(1))]
  fit <- do.call(ranger::ranger, args)

  pred <- stats::predict(
    fit,
    data = eval_df,
    num.threads = params$num.threads %||% 1L,
    verbose = FALSE
  )$predictions
  .drgpsm_prob_matrix(pred, tr_levels, "ranger GPS")
}

.drgpsm_fit_xgboost_outcome <- function(df, x_names, y, is_binary,
                                        params = NULL, seed = 12345) {
  .drgpsm_require_optional("xgboost", "outcome_model='xgboost'")

  x <- .drgpsm_model_matrix(df, x_names)
  defaults <- .drgpsm_xgboost_outcome_defaults(is_binary, seed)
  params <- utils::modifyList(defaults, params %||% list())

  fit <- .drgpsm_xgboost_train(
    x = x,
    y = y,
    params = params,
    objective = if (is_binary) "binary:logistic" else "reg:squarederror",
    eval_metric = if (is_binary) "logloss" else "rmse"
  )

  structure(
    list(model = fit, columns = colnames(x), is_binary = is_binary),
    class = "drgpsm_xgboost_outcome"
  )
}

.drgpsm_predict_xgboost_outcome <- function(fit, newdata, x_names) {
  x <- .drgpsm_model_matrix(newdata, x_names, columns = fit$columns)
  as.numeric(stats::predict(fit$model, newdata = xgboost::xgb.DMatrix(data = x)))
}

.drgpsm_ranger_outcome_defaults <- function(p, is_binary, seed) {
  list(
    num.trees = 1000L,
    mtry = if (is_binary) max(1L, floor(sqrt(p))) else max(1L, floor(p / 3)),
    min.node.size = if (is_binary) 10L else 5L,
    sample.fraction = 0.85,
    replace = FALSE,
    respect.unordered.factors = "order",
    oob.error = FALSE,
    num.threads = 1L,
    seed = seed,
    verbose = FALSE
  )
}

.drgpsm_fit_ranger_outcome <- function(df, x_names, outcome_var, is_binary,
                                       params = NULL, seed = 12345) {
  .drgpsm_require_optional("ranger", "outcome_model='ranger'")

  defaults <- .drgpsm_ranger_outcome_defaults(length(x_names), is_binary, seed)
  params <- utils::modifyList(defaults, params %||% list())
  params$probability <- NULL

  form <- stats::reformulate(x_names, response = outcome_var)
  args <- c(
    list(formula = form, data = df, probability = is_binary),
    params
  )
  args <- args[!vapply(args, is.null, logical(1))]
  fit <- do.call(ranger::ranger, args)

  structure(
    list(model = fit, is_binary = is_binary, num.threads = params$num.threads %||% 1L),
    class = "drgpsm_ranger_outcome"
  )
}

.drgpsm_predict_ranger_outcome <- function(fit, newdata) {
  pred <- stats::predict(
    fit$model,
    data = newdata,
    num.threads = fit$num.threads,
    verbose = FALSE
  )$predictions

  if (isTRUE(fit$is_binary)) {
    return(.drgpsm_prob_one(pred))
  }
  as.numeric(pred)
}
