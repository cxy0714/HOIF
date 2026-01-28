compute_HOIF <- function(X, A, Y, mu1, mu0, pi, m, sample_splitting = 0, K = 5) {
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
  compute_HOIF_a <- function(a, X_est, A_est, Y_est, R_a_est, s_a_est, r_a_est, Z_est, Omega_a) {
    n_est <- length(Y_est)

    # Generate all m-tuples of distinct indices
    # For small n, we enumerate all combinations
    indices <- 1:n_est

    # Generate all m-combinations without replacement
    if (m == 1) {
      tuples <- matrix(indices, ncol = 1)
    } else {
      tuples <- combn(indices, m)
      tuples <- t(tuples)  # Each row is one m-tuple
    }

    n_tuples <- nrow(tuples)

    # Sum over all m-tuples
    total_sum <- 0

    for (k in 1:n_tuples) {
      idx <- tuples[k, ]  # (i_1, i_2, ..., i_m)

      # Start with R^a_{i_1} * Z_{i_1}^T * Omega^a
      term <- R_a_est[idx[1]] * (t(Z_est[idx[1], , drop = FALSE]) %*% Omega_a)

      # Middle products for s = 2 to m-1
      if (m >= 3) {
        for (s in 2:(m-1)) {
          Q_a_is <- s_a_est[idx[s]] * Z_est[idx[s], , drop = FALSE] %*% t(Z_est[idx[s], , drop = FALSE])

          # Compute Sigma^a from estimation sample
          Sigma_a_est <- matrix(0, p, p)
          for (j in 1:n_est) {
            Sigma_a_est <- Sigma_a_est + s_a_est[j] * Z_est[j, , drop = FALSE] %*% t(Z_est[j, , drop = FALSE])
          }
          Sigma_a_est <- Sigma_a_est / n_est

          term <- term %*% (Q_a_is - Sigma_a_est) %*% Omega_a
        }
      }

      # Final term: Z_{i_m} * r^a_{i_m}
      term <- term %*% Z_est[idx[m], , drop = FALSE] * r_a_est[idx[m]]

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

    HOIF_1_m <- compute_HOIF_a(1, X, A, Y, R_1, s_1, r_1, Z, Omega_1)

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

    HOIF_0_m <- compute_HOIF_a(0, X, A, Y, R_0, s_0, r_0, Z, Omega_0)

  } else {
    # Sample splitting with K folds

    # Create fold indices
    folds <- cut(seq(1, n), breaks = K, labels = FALSE)
    folds <- sample(folds)  # Randomly assign to folds

    HOIF_1_m <- 0
    HOIF_0_m <- 0

    for (j in 1:K) {
      # Estimation sample (fold j)
      est_idx <- which(folds == j)
      # Training sample (all other folds)
      train_idx <- which(folds != j)

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

      HOIF_1_m <- HOIF_1_m + compute_HOIF_a(1, X_est, A_est, Y_est, R_1_est, s_1_est, r_1_est, Z_est, Omega_1_train)

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

      HOIF_0_m <- HOIF_0_m + compute_HOIF_a(0, X_est, A_est, Y_est, R_0_est, s_0_est, r_0_est, Z_est, Omega_0_train)
    }

    # Average over folds (optional, depending on your definition)
    # HOIF_1_m <- HOIF_1_m / K
    # HOIF_0_m <- HOIF_0_m / K
  }

  return(list(HOIF_1_m = HOIF_1_m, HOIF_0_m = HOIF_0_m))
}
