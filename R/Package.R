#' Package-level Documentation and Startup
#'
#' @author Xingyu Chen
#' @date 2026-01-23

# Package environment to store state
hoif_env <- new.env(parent = emptyenv())


#' Package startup function
#'
#' @param libname Library name
#' @param pkgname Package name
#'
#' @keywords internal
.onLoad <- function(libname, pkgname) {
  # Initialize cached status as NULL
  hoif_env$python_available <- NULL
  hoif_env$ustats_available <- NULL
}


#' HOIF: Higher Order Influence Function Estimators for ATE
#'
#' This package implements Higher Order Influence Function (HOIF) estimators
#' for Average Treatment Effect (ATE) estimation.
#'
#' @section Getting Started:
#' Before using this package, you need to set up the Python environment:
#'
#' \code{setup_hoif()}
#'
#' This will install Python (if needed) and the required u-stats package.
#'
#' @section Checking Setup:
#' To verify your environment is correctly configured:
#'
#' \code{check_hoif_setup()}
#'
#' @section Main Functions:
#' \itemize{
#'   \item \code{\link{hoif_ate}}: Main function to compute HOIF estimators
#'   \item \code{\link{ustat}}: Compute U-statistics (requires Python)
#'   \item \code{\link{setup_hoif}}: Setup Python environment
#'   \item \code{\link{check_hoif_setup}}: Check environment status
#' }
#'
#' @section Example Workflow:
#' \preformatted{
#' # 1. First-time setup (only once)
#' library(HOIF)
#' setup_hoif()
#'
#' # 2. Verify setup
#' check_hoif_setup()
#'
#' # 3. Use the package
#' results <- hoif_ate(X, A, Y, mu1, mu0, pi, ...)
#' print(results)
#' plot(results)
#' }
#'
#' @docType package
#' @name HOIF-package
#' @aliases HOIF
NULL
