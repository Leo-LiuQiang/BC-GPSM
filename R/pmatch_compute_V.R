#' Compute log-ratio PS index V from multinomial probabilities
#'
#' @param e_mat n x K matrix of probabilities (rows sum to 1)
#' @param ref_col Column index or column name of the reference treatment.
#'   Defaults to the last column.
#' @param eps small clipping value
#' @return n x (K-1) matrix of log-ratios relative to the reference column.
#' @export
#'
#' @examples
#' e_mat <- matrix(
#'   c(0.20, 0.50, 0.30,
#'     0.40, 0.40, 0.20),
#'   nrow = 2,
#'   byrow = TRUE
#' )
#' pmatch_compute_V(e_mat)
#' pmatch_compute_V(e_mat, ref_col = 1)
pmatch_compute_V <- function(e_mat, ref_col = ncol(e_mat), eps = 1e-12) {
  e_mat <- as.matrix(e_mat)
  if (ncol(e_mat) < 2L) stop("`e_mat` must have at least two probability columns.")
  if (is.character(ref_col)) {
    if (is.null(colnames(e_mat)) || !ref_col %in% colnames(e_mat)) {
      stop("Character `ref_col` must match a column name in `e_mat`.")
    }
    ref_col <- match(ref_col, colnames(e_mat))
  }
  if (!is.numeric(ref_col) || length(ref_col) != 1L ||
      ref_col < 1L || ref_col > ncol(e_mat)) {
    stop("`ref_col` must be a valid single column index or column name.")
  }
  e_mat <- pmax(e_mat, eps)
  e_mat <- e_mat / rowSums(e_mat)
  K <- ncol(e_mat)
  keep <- setdiff(seq_len(K), as.integer(ref_col))
  V <- log(e_mat[, keep, drop = FALSE] / e_mat[, ref_col])
  if (!is.null(colnames(e_mat))) {
    colnames(V) <- paste0(colnames(e_mat)[keep], "_vs_", colnames(e_mat)[ref_col])
  }
  V
}
