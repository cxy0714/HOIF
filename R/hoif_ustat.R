#' Compute a Higher-Order U-Statistic
#'
#' Computes a higher-order U-statistic given pre-computed kernel matrices
#' using an Einstein summation expression.
#'
#' @param tensors A list of numeric matrices of equal dimensions.
#' @param expression A character string specifying the Einstein summation.
#' @param backend Character string, either "numpy" or "torch".
#' @param average Logical; whether to return the averaged U-statistic.
#'
#' @return A numeric scalar.
#'
#' @examples
#' \dontrun{
#' H1 <- matrix(runif(100), 10, 10)
#' H2 <- matrix(runif(100), 10, 10)
#' hoif_ustat(list(H1, H2), "ab,bc->")
#' }
#'
#' @export
hoif_ustat <- function(tensors,
                       expression,
                       backend = c("numpy", "torch"),
                       average = TRUE) {

  backend <- match.arg(backend)

  if (!reticulate::py_available(initialize = FALSE)) {
    stop(
      "Python is required for hoif_ustat(). ",
      "Please install Python and the 'u-stats' package.",
      call. = FALSE
    )
  }

  if (backend == "torch" &&
      !reticulate::py_module_available("torch")) {
    warning(
      "Torch backend not available; falling back to numpy.",
      call. = FALSE
    )
    backend <- "numpy"
  }

  ustats <- reticulate::import("u_stats", delay_load = TRUE)
  ustats$set_backend(backend)

  np <- reticulate::import("numpy")

  tensors <- lapply(tensors, function(x) {
    if (!is.matrix(x)) {
      stop("All elements of 'tensors' must be matrices.", call. = FALSE)
    }
    np$array(x, dtype = "float32")
  })

  ustats$ustat(
    tensors = tensors,
    expression = expression,
    average = average
  )
}
