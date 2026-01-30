test_that("HOIF with sample splitting matches brute-force implementation", {

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

  m <- 7
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
    backend = "torch"
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
test_that("HOIF without sample splitting matches brute-force implementation", {

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

  m <- 7


  results <- hoif_ate(
    X, A, Y,
    mu1 = mu1,
    mu0 = mu0,
    pi = pi,
    transform_method = "none",
    m = m,
    sample_split = FALSE,
    backend = "torch"
  )

  HOIF_test <- HOIF:::compute_HOIF_sequence_test(
    X, A, Y, mu1, mu0, pi,
    m = m,
    sample_splitting = 0
  )

  expect_equal(results$HOIF1, HOIF_test$HOIF_1_m, tolerance = 1e-8)
  expect_equal(results$HOIF0, HOIF_test$HOIF_0_m, tolerance = 1e-8)
})
