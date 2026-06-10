# Pure R functionality: these tests do not require Python and therefore
# always run, including on CRAN.

test_that("transform_covariates returns correct dimensions", {
  set.seed(1)
  n <- 30
  p <- 2
  X <- matrix(rnorm(n * p), ncol = p)

  # method = "none": raw covariates, unchanged
  Z_none <- transform_covariates(X, method = "none")
  expect_equal(dim(Z_none), c(n, p))
  expect_equal(Z_none, X)

  # method = "splines": intercept + (basis_dim - 1) columns per covariate
  basis_dim <- 5
  Z_spl <- transform_covariates(X, method = "splines", basis_dim = basis_dim)
  expect_equal(nrow(Z_spl), n)
  expect_equal(ncol(Z_spl), 1 + p * (basis_dim - 1))
  expect_true(all(Z_spl[, 1] == 1))

  # method = "fourier": intercept + basis_dim columns per covariate
  Z_fou <- transform_covariates(X, method = "fourier", basis_dim = 4)
  expect_equal(nrow(Z_fou), n)
  expect_equal(ncol(Z_fou), 1 + p * 4)

  # odd basis_dim adds one extra cosine column
  Z_fou_odd <- transform_covariates(X, method = "fourier", basis_dim = 5)
  expect_equal(ncol(Z_fou_odd), 1 + p * 5)

  # informative errors
  expect_error(transform_covariates(X, method = "splines"), "basis_dim")
  expect_error(transform_covariates(X, method = "unknown", basis_dim = 4),
               "Method must be")
})

test_that("compute_residuals returns expected values", {
  A <- c(1, 0, 1, 0)
  Y <- c(2, 1, 3, 0)
  mu1 <- rep(1, 4)
  mu0 <- rep(0.5, 4)
  pi <- rep(0.5, 4)

  res <- compute_residuals(A, Y, mu1, mu0, pi)

  expect_equal(res$R1, Y - mu1)
  expect_equal(res$R0, Y - mu0)
  expect_equal(res$r1, 1 - A / pi)
  expect_equal(res$r0, 1 - (1 - A) / (1 - pi))
})

test_that("build_Ej constructs the chain expression", {
  E3 <- build_Ej(3)
  expect_equal(E3, list(1, c(1, 2), c(2, 3), 3))
  expect_error(build_Ej(1), "j must be >= 2")
})

test_that("pure R HOIF (no sample splitting) matches brute-force implementation", {

  set.seed(42)

  n <- 7
  p <- 1

  X <- matrix(rnorm(n * p), ncol = p)
  A <- rbinom(n, 1, 0.5)

  beta <- runif(p)
  beta <- beta / sqrt(as.numeric(crossprod(beta)))

  Y <- as.numeric(A + X %*% beta + rnorm(n, 0, 0.1))
  mu1 <- as.numeric(1 + X %*% beta + rnorm(n, 0, 0.1))
  mu0 <- as.numeric(0 + X %*% beta + rnorm(n, 0, 0.1))
  pi <- rep(0.5, n)

  m <- 6

  results <- hoif_ate(
    X, A, Y,
    mu1 = mu1,
    mu0 = mu0,
    pi = pi,
    transform_method = "none",
    m = m,
    sample_split = FALSE,
    pure_R_code = TRUE
  )

  HOIF_test <- HOIF:::compute_HOIF_sequence_test(
    X, A, Y, mu1, mu0, pi,
    m = m,
    sample_splitting = 0
  )

  expect_equal(results$HOIF1, HOIF_test$HOIF_1_m, tolerance = 1e-8)
  expect_equal(results$HOIF0, HOIF_test$HOIF_0_m, tolerance = 1e-8)
})

test_that("pure R HOIF (with sample splitting) matches brute-force implementation", {

  skip_on_cran() # brute-force enumeration is slow

  set.seed(42)

  n <- 14
  p <- 1

  X <- matrix(rnorm(n * p), ncol = p)
  A <- rbinom(n, 1, 0.5)

  beta <- runif(p)
  beta <- beta / sqrt(as.numeric(crossprod(beta)))

  Y <- as.numeric(A + X %*% beta + rnorm(n, 0, 0.1))
  mu1 <- as.numeric(1 + X %*% beta + rnorm(n, 0, 0.1))
  mu0 <- as.numeric(0 + X %*% beta + rnorm(n, 0, 0.1))
  pi <- rep(0.5, n)

  m <- 6
  n_folds <- 2
  seed <- 42

  results_split <- hoif_ate(
    X, A, Y,
    mu1 = mu1,
    mu0 = mu0,
    pi = pi,
    transform_method = "none",
    m = m,
    sample_split = TRUE,
    n_folds = n_folds,
    seed = seed,
    pure_R_code = TRUE
  )

  HOIF_test_split <- HOIF:::compute_HOIF_sequence_test(
    X, A, Y, mu1, mu0, pi,
    m = m,
    sample_splitting = 1,
    n_folds = n_folds,
    seed = seed
  )

  expect_equal(results_split$HOIF1, HOIF_test_split$HOIF_1_m, tolerance = 1e-8)
  expect_equal(results_split$HOIF0, HOIF_test_split$HOIF_0_m, tolerance = 1e-8)
})

test_that("pure R mode enforces the order limit", {
  set.seed(7)
  n <- 30
  X <- matrix(rnorm(n), ncol = 1)
  A <- rbinom(n, 1, 0.5)
  Y <- rnorm(n)
  mu1 <- rep(0, n)
  mu0 <- rep(0, n)
  pi <- rep(0.5, n)

  expect_error(
    hoif_ate(X, A, Y, mu1 = mu1, mu0 = mu0, pi = pi,
             m = 7, pure_R_code = TRUE),
    "maximum order supported is 6"
  )
  expect_warning(
    hoif_ate(X, A, Y, mu1 = mu1, mu0 = mu0, pi = pi,
             m = 3, pure_R_code = TRUE),
    "up to 6"
  )
})

test_that("print and plot methods work", {
  set.seed(11)
  n <- 60
  X <- matrix(rnorm(n), ncol = 1)
  A <- rbinom(n, 1, 0.5)
  Y <- as.numeric(A + X[, 1] + rnorm(n, 0, 0.1))
  mu1 <- 1 + X[, 1]
  mu0 <- X[, 1]
  pi <- rep(0.5, n)

  fit <- hoif_ate(X, A, Y, mu1 = mu1, mu0 = mu0, pi = pi,
                  m = 6, pure_R_code = TRUE)

  expect_s3_class(fit, "hoif_ate")
  out <- capture.output(expect_invisible(print(fit)))
  expect_true(any(grepl("HOIF Estimators", out)))

  pdf(NULL)
  on.exit(dev.off(), add = TRUE)
  expect_silent(plot(fit))
})
