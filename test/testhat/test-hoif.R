#' Unit Tests for HOIF Estimators for ATE
#'
#' This file contains comprehensive unit tests for the HOIF package.
#' Tests use the testthat framework.
#' Tests that require Python are skipped if Python is not available.
#'
#' @author Xingyu Chen
#' @date 2026-01-23

library(testthat)

# Helper function to skip tests if Python not available
skip_if_no_python <- function() {
  if (!check_python_env()) {
    skip("Python or u-stats not available")
  }
}

# ==============================================================================
# Test 1: Expression conversion
# ==============================================================================

test_that("expr_list_to_einstein converts nested lists correctly", {
  # Test case 1: Simple case
  expr1 <- list(c(1,2), c(2,3))
  result1 <- expr_list_to_einstein(expr1)
  expect_equal(result1, "ab,bc->")

  # Test case 2: Longer chain
  expr2 <- list(c(1,2), c(2,3), c(3,4))
  result2 <- expr_list_to_einstein(expr2)
  expect_equal(result2, "ab,bc,cd->")

  # Test case 3: Non-sequential indices
  expr3 <- list(c(1,3), c(3,5), c(5,7))
  result3 <- expr_list_to_einstein(expr3)
  expect_true(grepl("->$", result3))  # Should end with ->
  expect_equal(length(strsplit(result3, ",")[[1]]), 3)  # Should have 3 terms
})


# ==============================================================================
# Test 2: Covariate transformation
# ==============================================================================

test_that("transform_covariates works with splines", {
  set.seed(123)
  n <- 100
  p <- 3
  X <- matrix(rnorm(n * p), ncol = p)

  # Test splines transformation
  Z <- transform_covariates(X, method = "splines", k = 5, degree = 3)

  expect_true(is.matrix(Z))
  expect_equal(nrow(Z), n)
  expect_true(ncol(Z) > p)  # Should expand dimension
  expect_equal(Z[, 1], rep(1, n))  # First column should be intercept
})

test_that("transform_covariates works with fourier", {
  set.seed(123)
  n <- 100
  p <- 2
  X <- matrix(rnorm(n * p), ncol = p)

  # Test Fourier transformation
  Z <- transform_covariates(X, method = "fourier", k = 6)

  expect_true(is.matrix(Z))
  expect_equal(nrow(Z), n)
  expect_equal(Z[, 1], rep(1, n))  # First column should be intercept
})

test_that("transform_covariates scales input correctly", {
  set.seed(123)
  X <- matrix(c(1, 2, 3, 10, 20, 30), ncol = 2)

  Z <- transform_covariates(X, method = "splines", k = 3)

  # All transformations should be based on [0,1] scaled input
  # Check that extreme values don't cause issues
  expect_true(all(is.finite(Z)))
})


# ==============================================================================
# Test 3: Residual computation
# ==============================================================================

test_that("compute_residuals returns correct structure", {
  set.seed(123)
  n <- 100
  A <- rbinom(n, 1, 0.5)
  Y <- rnorm(n)
  mu1 <- rnorm(n)
  mu0 <- rnorm(n)
  pi <- runif(n, 0.2, 0.8)

  residuals <- compute_residuals(A, Y, mu1, mu0, pi)

  expect_type(residuals, "list")
  expect_named(residuals, c("R1", "r1", "R0", "r0"))
  expect_equal(length(residuals$R1), n)
  expect_equal(length(residuals$r1), n)
  expect_equal(length(residuals$R0), n)
  expect_equal(length(residuals$r0), n)
})

test_that("compute_residuals calculates correctly", {
  n <- 10
  A <- c(1, 0, 1, 0, 1, 0, 1, 0, 1, 0)
  Y <- 1:10
  mu1 <- rep(5, n)
  mu0 <- rep(3, n)
  pi <- rep(0.5, n)

  residuals <- compute_residuals(A, Y, mu1, mu0, pi)

  # Check R1: should be A * (Y - mu1)
  expect_equal(residuals$R1[1], 1 * (1 - 5))  # A=1, Y=1, mu1=5
  expect_equal(residuals$R1[2], 0)  # A=0

  # Check r1: should be 1 - A/pi
  expect_equal(residuals$r1[1], 1 - 1/0.5)  # A=1, pi=0.5
  expect_equal(residuals$r1[2], 1 - 0/0.5)  # A=0

  # Check R0: should be (1-A) * (Y - mu0)
  expect_equal(residuals$R0[1], 0)  # A=1
  expect_equal(residuals$R0[2], 1 * (2 - 3))  # A=0, Y=2, mu0=3
})


# ==============================================================================
# Test 4: Gram matrix inversion
# ==============================================================================

test_that("compute_gram_inverse works with direct method", {
  set.seed(123)
  n <- 50
  k <- 10
  Z <- matrix(rnorm(n * k), ncol = k)
  A <- rbinom(n, 1, 0.5)

  Omega <- compute_gram_inverse(Z, A, method = "direct")

  expect_type(Omega, "list")
  expect_named(Omega, c("Omega1", "Omega0"))
  expect_true(is.matrix(Omega$Omega1))
  expect_true(is.matrix(Omega$Omega0))
  expect_equal(dim(Omega$Omega1), c(k, k))
  expect_equal(dim(Omega$Omega0), c(k, k))
})

test_that("compute_gram_inverse handles singular matrices", {
  set.seed(123)
  n <- 10
  k <- 20  # k > n will likely cause singularity
  Z <- matrix(rnorm(n * k), ncol = k)
  A <- rbinom(n, 1, 0.5)

  # Should not error, should use fallback
  expect_warning(
    Omega <- compute_gram_inverse(Z, A, method = "direct"),
    "not positive definite"
  )

  expect_true(is.matrix(Omega$Omega1))
  expect_true(is.matrix(Omega$Omega0))
})


# ==============================================================================
# Test 5: Basis matrix computation
# ==============================================================================

test_that("compute_basis_matrix returns correct dimensions", {
  set.seed(123)
  n <- 50
  k <- 10
  Z <- matrix(rnorm(n * k), ncol = k)
  Omega1 <- diag(k)
  Omega0 <- diag(k)

  B_matrices <- compute_basis_matrix(Z, Omega1, Omega0)

  expect_type(B_matrices, "list")
  expect_named(B_matrices, c("B1", "B0"))
  expect_equal(dim(B_matrices$B1), c(n, n))
  expect_equal(dim(B_matrices$B0), c(n, n))
})

test_that("compute_basis_matrix is symmetric", {
  set.seed(123)
  n <- 30
  k <- 5
  Z <- matrix(rnorm(n * k), ncol = k)
  Omega1 <- crossprod(matrix(rnorm(k * k), k, k))  # Positive definite
  Omega0 <- crossprod(matrix(rnorm(k * k), k, k))

  B_matrices <- compute_basis_matrix(Z, Omega1, Omega0)

  # Projection matrices should be symmetric
  expect_true(all(abs(B_matrices$B1 - t(B_matrices$B1)) < 1e-10))
  expect_true(all(abs(B_matrices$B0 - t(B_matrices$B0)) < 1e-10))
})


# ==============================================================================
# Test 6: HOIF estimators (integration test with mock ustat)
# ==============================================================================

test_that("compute_hoif_estimators returns correct structure", {
  skip_if_no_python()

  set.seed(123)
  n <- 50

  # Create simple test data
  residuals <- list(
    R1 = rnorm(n),
    r1 = rnorm(n),
    R0 = rnorm(n),
    r0 = rnorm(n)
  )

  B_matrices <- list(
    B1 = diag(n),  # Identity for simplicity
    B0 = diag(n)
  )

  m <- 4
  results <- compute_hoif_estimators(residuals, B_matrices, m = m,
                                     backend = "numpy")

  expect_type(results, "list")
  expect_named(results, c("ATE", "HOIF1", "HOIF0", "IIFF1", "IIFF0", "orders"))
  expect_equal(length(results$ATE), m - 1)
  expect_equal(results$orders, 2:m)
})


# ==============================================================================
# Test 7: Main function hoif_ate
# ==============================================================================

test_that("hoif_ate works without sample splitting", {
  skip_if_no_python()

  set.seed(123)
  n <- 100
  p <- 3

  # Generate simple data
  X <- matrix(rnorm(n * p), ncol = p)
  A <- rbinom(n, 1, 0.5)
  Y <- A + X[, 1] + rnorm(n, sd = 0.5)

  # Simple nuisance estimates
  mu1 <- rep(mean(Y[A == 1]), n)
  mu0 <- rep(mean(Y[A == 0]), n)
  pi <- rep(0.5, n)

  # Run HOIF
  results <- hoif_ate(
    X = X, A = A, Y = Y, mu1 = mu1, mu0 = mu0, pi = pi,
    transform_method = "splines", k = 5, m = 3,
    sample_split = FALSE, backend = "numpy"
  )

  expect_s3_class(results, "hoif_ate")
  expect_true(is.numeric(results$ATE))
  expect_equal(length(results$ATE), 2)  # m=3 means orders 2,3
})

test_that("hoif_ate works with sample splitting", {
  skip_if_no_python()

  set.seed(456)
  n <- 100
  p <- 2

  X <- matrix(rnorm(n * p), ncol = p)
  A <- rbinom(n, 1, 0.5)
  Y <- A + X[, 1] + rnorm(n, sd = 0.5)

  mu1 <- rep(mean(Y[A == 1]), n)
  mu0 <- rep(mean(Y[A == 0]), n)
  pi <- rep(0.5, n)

  # Run with sample splitting
  results <- hoif_ate(
    X = X, A = A, Y = Y, mu1 = mu1, mu0 = mu0, pi = pi,
    transform_method = "splines", k = 5, m = 3,
    sample_split = TRUE, K = 3, seed = 789,
    backend = "numpy"
  )

  expect_s3_class(results, "hoif_ate")
  expect_true(is.numeric(results$ATE))

  # Test reproducibility with seed
  results2 <- hoif_ate(
    X = X, A = A, Y = Y, mu1 = mu1, mu0 = mu0, pi = pi,
    transform_method = "splines", k = 5, m = 3,
    sample_split = TRUE, K = 3, seed = 789,
    backend = "numpy"
  )

  expect_equal(results$ATE, results2$ATE)
})


# ==============================================================================
# Test 8: Print and plot methods
# ==============================================================================

test_that("print.hoif_ate works", {
  # Create mock results
  mock_results <- list(
    ATE = c(0.5, 0.52),
    HOIF1 = c(1.0, 1.05),
    HOIF0 = c(0.5, 0.53),
    IIFF1 = c(0.3, 0.05),
    IIFF0 = c(0.2, 0.03),
    orders = 2:3,
    convergence_data = data.frame(order = 2:3, ATE = c(0.5, 0.52))
  )
  class(mock_results) <- "hoif_ate"

  # Should not error
  expect_output(print(mock_results), "HOIF Estimators")
  expect_output(print(mock_results), "ATE")
})

test_that("plot.hoif_ate works", {
  mock_results <- list(
    ATE = c(0.5, 0.52, 0.53),
    orders = 2:4,
    convergence_data = data.frame(order = 2:4, ATE = c(0.5, 0.52, 0.53))
  )
  class(mock_results) <- "hoif_ate"

  # Should not error
  expect_silent(plot(mock_results))
})


# ==============================================================================
# Test 9: Edge cases and error handling
# ==============================================================================

test_that("hoif_ate handles edge cases", {
  skip_if_no_python()

  n <- 50

  # Test with single covariate
  X <- matrix(rnorm(n), ncol = 1)
  A <- rbinom(n, 1, 0.5)
  Y <- A + rnorm(n)
  mu1 <- rep(0, n)
  mu0 <- rep(0, n)
  pi <- rep(0.5, n)

  expect_silent(
    results <- hoif_ate(X, A, Y, mu1, mu0, pi, k = 3, m = 2,
                        sample_split = FALSE, backend = "numpy")
  )
})

test_that("transform_covariates errors on invalid method", {
  X <- matrix(rnorm(20), ncol = 2)

  expect_error(
    transform_covariates(X, method = "invalid"),
    "Method must be 'splines' or 'fourier'"
  )
})

test_that("compute_gram_inverse errors on invalid method", {
  Z <- matrix(rnorm(50), ncol = 5)
  A <- rbinom(10, 1, 0.5)

  expect_error(
    compute_gram_inverse(Z, A, method = "invalid"),
    "Method must be"
  )
})


# ==============================================================================
# Test 10: Numerical accuracy tests
# ==============================================================================

test_that("Residuals sum to correct values under null", {
  set.seed(999)
  n <- 200

  # Under null: no treatment effect, correctly specified
  A <- rbinom(n, 1, 0.5)
  X <- rnorm(n)
  Y <- X + rnorm(n, sd = 0.1)  # No treatment effect

  # Perfect nuisance estimates
  mu1 <- X
  mu0 <- X
  pi <- rep(0.5, n)

  residuals <- compute_residuals(A, Y, mu1, mu0, pi)

  # Under correct specification, residuals should be small
  expect_true(abs(mean(residuals$R1)) < 0.1)
  expect_true(abs(mean(residuals$R0)) < 0.1)
})


# ==============================================================================
# Run all tests
# ==============================================================================

# To run these tests, use:
# testthat::test_file("tests/testthat/test-hoif.R")
# Or in the package: devtools::test()
