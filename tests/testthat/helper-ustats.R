# Helper: check whether the Python backend (u_stats) is available without
# triggering a download on machines where it is not set up.
has_ustats <- function() {
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    return(FALSE)
  }
  isTRUE(tryCatch(
    reticulate::py_available(initialize = TRUE) &&
      reticulate::py_module_available("u_stats"),
    error = function(e) FALSE
  ))
}

skip_if_no_ustats <- function() {
  # Never initialize Python (or trigger automatic dependency downloads)
  # on CRAN machines.
  testthat::skip_on_cran()
  if (!has_ustats()) {
    testthat::skip("Python/u_stats not available")
  }
}
