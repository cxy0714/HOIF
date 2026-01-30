
# Helper function to generate all m-permutations (ordered tuples) of distinct indices
generate_permutations <- function(indices, m) {
  n_idx <- length(indices)
  if (m == 1) {
    return(matrix(indices, ncol = 1))
  }

  # Use recursive approach or gtools::permutations
  # For simplicity, here's a manual approach:
  # First get all m-combinations, then permute each

  if (m > n_idx) {
    stop("m cannot be larger than the number of indices")
  }

  # Get all m-combinations
  combs <- combn(indices, m, simplify = FALSE)

  # For each combination, generate all permutations
  all_perms <- list()
  for (comb in combs) {
    # Generate all permutations of this combination
    perms <- permn(comb)  # This will use a helper function
    all_perms <- c(all_perms, perms)
  }

  # Convert list to matrix
  result <- do.call(rbind, lapply(all_perms, function(x) matrix(x, nrow = 1)))
  return(result)
}

# Helper function to generate all permutations of a vector
permn <- function(x) {
  if (length(x) == 1) {
    return(list(x))
  }
  result <- list()
  for (i in seq_along(x)) {
    rest <- x[-i]
    for (p in permn(rest)) {
      result <- c(result, list(c(x[i], p)))
    }
  }
  return(result)
}

compute_HOIF_sequence <- function(X, A, Y, mu1, mu0, pi, m, sample_splitting = 0, K = 2, seed) {
  HOIF_1_m <- numeric(m-1)
  HOIF_0_m <- numeric(m-1)
  IIFF_1_m <- numeric(m-1)
  IIFF_0_m <- numeric(m-1)
  for (i in 2:m)
  {
    results <- compute_HOIF(X, A, Y, mu1, mu0, pi, i, sample_splitting, K,seed )
    IIFF_1_m[i -1] <- results$IIFF_1_m
    IIFF_0_m[i -1] <- results$IIFF_0_m
  }
  for ( i in 2:m) {
    for ( j in 2 : i){

    HOIF_1_m[i -1] <- HOIF_1_m[i -1] + IIFF_1_m[j - 1]
    HOIF_0_m[i -1] <- HOIF_0_m[i -1] + IIFF_0_m[j - 1]
    }
  }
  return(
    list(
      HOIF_1_m = HOIF_1_m,
      HOIF_0_m = HOIF_0_m,
      IIFF_1_m = IIFF_1_m,
      IIFF_0_m = IIFF_0_m
    )
  )
}
compute_HOIF <- function(X, A, Y, mu1, mu0, pi, m, sample_splitting = 0, K = 2, seed) {
  # X: n x p covariate matrix
  # A: n-vector of binary treatment (0 or 1)
  # Y: n-vector of outcomes
  # mu1: n-vector of estimated mu(1, X_i)
  # mu0: n-vector of estimated mu(0, X_i)
  # pi: n-vector of estimated propensity scores
  # m: order of the HOIF statistic
  # sample_splitting: 0 for no splitting, 1 for K-fold splitting
  # K: number of folds (only used when sample_splitting = 1)

  n <- length(Y)
  p <- ncol(X)



  # Helper function to compute HOIF for a given a (treatment arm)
  compute_HOIF_a <- function( X_est, A_est, Y_est, R_a_est, s_a_est, r_a_est, Z_est, Omega_a, Sigma_a) {
    n_est <- length(Y_est)
    p <- ncol(Z_est)

    # Generate all m-tuples of distinct indices (ordered)
    indices <- 1:n_est
    tuples <- generate_permutations(indices, m)
    n_tuples <- nrow(tuples)

    # Sum over all m-tuples
    total_sum <- 0

    # Compute Sigma^a from estimation sample

    for (k in 1:n_tuples) {
      idx <- tuples[k, ]  # (i_1, i_2, ..., i_m)

      # Start with R^a_{i_1} * Z_{i_1}^T * Omega^a
      term <- R_a_est[idx[1]] * (t(Z_est[idx[1], , drop = FALSE]) %*% Omega_a)

      # Middle products for s = 2 to m-1
      if (m >= 3) {
        for (s in 2:(m-1)) {
          Q_a_is <- s_a_est[idx[s]] * Z_est[idx[s], , drop = FALSE] %*% t(Z_est[idx[s], , drop = FALSE])



          term <- term %*% (Q_a_is - Sigma_a) %*% Omega_a
          # term <- term %*% (Q_a_is ) %*% Omega_a
        }
      }

      # Final term: Z_{i_m} * r^a_{i_m}
      term <- term %*% Z_est[idx[m], , drop = FALSE] * r_a_est[idx[m]] *  s_a_est[idx[m]]

      total_sum <- total_sum + term[1, 1]
    }

    # Multiply by coefficient
    coef <- (-1)^m * factorial(n_est - m) / factorial(n_est)

    return(coef * total_sum)
  }

  # Main computation
  if (sample_splitting == 0) {
    # No sample splitting

    # Compute for a = 1
    R_1 <- Y - mu1
    s_1 <- A
    r_1 <- 1 - s_1 / pi
    Z <- X

    # Compute Sigma^1 and Omega^1
    Sigma_1 <- matrix(0, p, p)
    for (i in 1:n) {
      Sigma_1 <- Sigma_1 + s_1[i] * Z[i, , drop = FALSE] %*% t(Z[i, , drop = FALSE])
    }
    Sigma_1 <- Sigma_1 / n
    Omega_1 <- solve(Sigma_1)

    IIFF_1_m <- compute_HOIF_a( X, A, Y, R_1, s_1, r_1, Z, Omega_1, Sigma_1)

    # Compute for a = 0
    R_0 <- Y - mu0
    s_0 <- 1 - A
    r_0 <- 1 - s_0 / (1 - pi)

    # Compute Sigma^0 and Omega^0
    Sigma_0 <- matrix(0, p, p)
    for (i in 1:n) {
      Sigma_0 <- Sigma_0 + s_0[i] * Z[i, , drop = FALSE] %*% t(Z[i, , drop = FALSE])
    }
    Sigma_0 <- Sigma_0 / n
    Omega_0 <- solve(Sigma_0)

    IIFF_0_m <- compute_HOIF_a( X, A, Y, R_0, s_0, r_0, Z, Omega_0, Sigma_0)

  } else {
    # Sample splitting with K folds
    set.seed(seed)
    # Create fold indices
    fold_indices <- sample(rep(1:n_folds, length.out = n))

    IIFF_1_m <- numeric(K)
    IIFF_0_m <- numeric(K)

    for (j in 1:K) {
      # Estimation sample (fold j)
      est_idx <- which(fold_indices == j)
      # Training sample (all other folds)
      train_idx <- which(fold_indices != j)

      # Training data
      X_train <- X[train_idx, , drop = FALSE]
      A_train <- A[train_idx]
      Y_train <- Y[train_idx]
      mu1_train <- mu1[train_idx]
      mu0_train <- mu0[train_idx]
      pi_train <- pi[train_idx]

      # Estimation data
      X_est <- X[est_idx, , drop = FALSE]
      A_est <- A[est_idx]
      Y_est <- Y[est_idx]
      mu1_est <- mu1[est_idx]
      mu0_est <- mu0[est_idx]
      pi_est <- pi[est_idx]

      Z_train <- X_train
      Z_est <- X_est

      # For a = 1
      s_1_train <- A_train
      s_1_est <- A_est
      R_1_est <- Y_est - mu1_est
      r_1_est <- 1 - s_1_est / pi_est

      # Compute Omega^1 from training data
      Sigma_1_train <- matrix(0, p, p)
      for (i in 1:length(train_idx)) {
        Sigma_1_train <- Sigma_1_train + s_1_train[i] * Z_train[i, , drop = FALSE] %*% t(Z_train[i, , drop = FALSE])
      }
      Sigma_1_train <- Sigma_1_train / length(train_idx)
      Omega_1_train <- solve(Sigma_1_train)

      IIFF_1_m[j] <-  compute_HOIF_a( X_est, A_est, Y_est, R_1_est, s_1_est, r_1_est, Z_est, Omega_1_train, Sigma_1_train)

      # For a = 0
      s_0_train <- 1 - A_train
      s_0_est <- 1 - A_est
      R_0_est <- Y_est - mu0_est
      r_0_est <- 1 - s_0_est / (1 - pi_est)

      # Compute Omega^0 from training data
      Sigma_0_train <- matrix(0, p, p)
      for (i in 1:length(train_idx)) {
        Sigma_0_train <- Sigma_0_train + s_0_train[i] * Z_train[i, , drop = FALSE] %*% t(Z_train[i, , drop = FALSE])
      }
      Sigma_0_train <- Sigma_0_train / length(train_idx)
      Omega_0_train <- solve(Sigma_0_train)

      IIFF_0_m[j] <-  compute_HOIF_a( X_est, A_est, Y_est, R_0_est, s_0_est, r_0_est, Z_est, Omega_0_train, Sigma_0_train)
    }
    IIFF_0_m <- mean(IIFF_0_m)
    IIFF_1_m <- mean(IIFF_1_m)
  }

  return(list(IIFF_1_m = IIFF_1_m, IIFF_0_m = IIFF_0_m))
}

devtools::load_all()

library(ustats)
library(HOIF)

# ==============================================================================
# Test 1: Environment Setup
# ==============================================================================
cat("\n========================================\n")
cat("Test 1: Check Environment of ustats \n")
cat("========================================\n")

library(reticulate)
py_config()

# setup_ustats()

check_ustats_setup()


set.seed(123)

n <- 10 # sample size
p <- 1 # number of covariates

# Covariates
X <- matrix(rnorm(n * p), ncol = p)

# Binary treatment assignment (randomized)
A <- rbinom(n, 1, 0.5)
# A <- rep(1,n)
# True coefficient vector (normalized to unit length)
beta <- runif(p)
beta <- beta / sqrt(as.numeric(crossprod(beta)))

# Outcome: partially linear model
Y <- as.numeric(A + X %*% beta + rnorm(n, 0, 0.1))


# ------------------------------------------------------------------------------
# Step 2: Fit nuisance outcome regressions via Lasso
# ------------------------------------------------------------------------------

mu1 <- as.numeric(1 + X %*% beta + rnorm(n, 0, 0.1))
mu0 <- as.numeric(0 + X %*% beta + rnorm(n, 0, 0.1))

# Known propensity score
pi <- rep(0.5, n)


# ------------------------------------------------------------------------------
# Step 3: Configure HOIF estimator with sample splitting
# ------------------------------------------------------------------------------
# m               : Order of the HOIF expansion (see ?hoif_ate for details)
# n_folds         : Number of folds used for cross-fitting
#                   n_folds = 2 corresponds to the eHOIF setting
# transform_method: Feature transformation applied to X before estimation.
#                   "none" means raw covariates are used (see ?hoif_ate)


m <- 5
n_folds <- 2
seed <- 123

# ------------------------------------------------------------------------------
# Step 4: Run HOIF
# ------------------------------------------------------------------------------


cat("\nResults (without sample splitting):\n")


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

cat("\nResults (with sample splitting):\n")
print(results_split)


HOIF_test_split <- compute_HOIF_sequence(X,A,Y,mu1,mu0, pi, m = m, sample_splitting = 1, K = n_folds,seed = seed)

all.equal(results_split$HOIF1, HOIF_test_split$HOIF_1_m)
results_split$HOIF1
HOIF_test_split$HOIF_1_m
all.equal(results_split$HOIF0, HOIF_test_split$HOIF_0_m)
results_split$HOIF0
HOIF_test_split$HOIF_0_m

results <- hoif_ate(
  X, A, Y,
  mu1 = mu1,
  mu0 = mu0,
  pi = pi,
  transform_method = "none",
  m = m,
  sample_split = FALSE,
  n_folds = n_folds,
  seed = seed,
  backend = "torch"
)

HOIF_test <- compute_HOIF_sequence(X,A,Y,mu1,mu0, pi, m = m, sample_splitting = 0, K = n_folds,seed = seed)


all.equal(results$HOIF1, HOIF_test$HOIF_1_m )
results$HOIF1
HOIF_test$HOIF_1_m

all.equal(results$HOIF0, HOIF_test$HOIF_0_m )
results$HOIF0
HOIF_test$HOIF_0_m
