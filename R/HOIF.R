#' HOIF Estimators for Average Treatment Effect
#'
#' This file implements the Higher Order Influence Function (HOIF) estimators
#' for Average Treatment Effect (ATE) estimation with nuisance functions.
#'
#' @author Xingyu Chen



#---------------------------core function----------------------------------
#' Transform covariates to basis functions
#'
#' @param X Matrix of covariates (n x p)
#' @param method Character: "splines", "fourier", or "none"
#' @param basis_dim Integer: dimension of basis expansion (ignored if method = "none")
#' @param degree Integer: degree for B-splines (default 3; ignored if method != "splines")
#' @param period Numeric: period for Fourier basis (default 1; ignored if method != "fourier")
#'
#' @return Matrix Z (n x (k * p + 1) or n x (p + 1)) of transformed covariates with intercept
#' @export
transform_covariates <- function(X, method = "splines", basis_dim, degree = 3, period = 1) {
  if (is.vector(X)) X <- matrix(X, ncol = 1)
  n <- nrow(X)
  p <- ncol(X)



  if (method == "none") {
    # Return [intercept | X] as-is
    return(X)
  }

  # Initialize with intercept
  Z_intercept <- matrix(1, nrow = n, ncol = 1)
  # Scale X to [0, 1] for each dimension (only needed for splines/fourier)
  X_scaled <- apply(X, 2, function(x) {
    (x - min(x)) / (max(x) - min(x) + 1e-10)
  })

  if (method == "splines") {
    Z_basis <- matrix(nrow = n, ncol = 0)
    for (j in 1:p) {
      knots <- seq(0, 1, length.out = k - degree + 1)
      basis_j <- splines::bs(
        X_scaled[, j],
        knots = knots[-c(1, length(knots))],
        degree = degree,
        intercept = FALSE
      )
      Z_basis <- cbind(Z_basis, basis_j)
    }
    return(cbind(Z_intercept, Z_basis))

  } else if (method == "fourier") {
    Z_basis <- matrix(nrow = n, ncol = 0)
    for (j in 1:p) {
      freq <- 1:floor(k / 2)
      for (f in freq) {
        Z_basis <- cbind(
          Z_basis,
          cos(2 * pi * f * X_scaled[, j] / period),
          sin(2 * pi * f * X_scaled[, j] / period)
        )
      }
      # Handle odd k: add one extra cosine if k is odd
      if (k %% 2 == 1) {
        f_extra <- (k + 1) / 2
        Z_basis <- cbind(Z_basis, cos(2 * pi * f_extra * X_scaled[, j] / period))
      }
    }
    return(cbind(Z_intercept, Z_basis))

  } else {
    stop("Method must be 'splines', 'fourier', or 'none'")
  }
}

#' Compute residuals for both treatment groups
#'
#' @param A Treatment vector (0 or 1)
#' @param Y Outcome vector
#' @param mu1 Predicted outcomes under treatment (mu(1, X))
#' @param mu0 Predicted outcomes under control (mu(0, X))
#' @param pi Propensity scores
#'
#' @return List with R1, r1, R0, r0
#' @export
compute_residuals <- function(A, Y, mu1, mu0, pi) {
  # For a = 1
  R1 <- Y - mu1
  r1 <- 1 - A / pi

  # For a = 0
  R0 <- Y - mu0
  r0 <- 1 - (1 - A) / (1 - pi)

  return(list(R1 = R1, r1 = r1, R0 = R0, r0 = r0))
}


#' Compute inverse of weighted Gram matrix
#'
#' @param Z Basis matrix (n x k)
#' @param A Treatment vector
#' @param method Character: "direct", "nlshrink", or "corpcor"
#'
#' @return List with Omega1 and Omega0 (inverse Gram matrices)
#' @export
compute_gram_inverse <- function(Z, A, method = "direct") {
  n <- nrow(Z)
  k <- ncol(Z)

  # Compute weights
  s1 <- A  # s_i^1 = A_i
  s0 <- 1 - A  # s_i^0 = 1 - A_i

  # Weighted Gram matrices using eigenMapMatMult if available
  if (requireNamespace("SMUT", quietly = TRUE)) {
    # Use faster matrix multiplication
    Z_weighted1 <- Z * sqrt(s1)
    Z_weighted0 <- Z * sqrt(s0)
    G1 <- SMUT::eigenMapMatMult(t(Z_weighted1), Z_weighted1) / n
    G0 <- SMUT::eigenMapMatMult(t(Z_weighted0), Z_weighted0) / n
  } else {
    # Fallback to standard crossprod
    G1 <- crossprod(Z * sqrt(s1)) / n
    G0 <- crossprod(Z * sqrt(s0)) / n
  }

  # Compute inverses
  if (method == "direct") {
    Omega1 <- tryCatch({
      chol2inv(chol(G1))
    }, error = function(e) {
      warning("G1 not positive definite, using Moore-Penrose inverse")
      MASS::ginv(G1)
    })

    Omega0 <- tryCatch({
      chol2inv(chol(G0))
    }, error = function(e) {
      warning("G0 not positive definite, using Moore-Penrose inverse")
      MASS::ginv(G0)
    })
  } else if (method == "nlshrink") {
    if (!requireNamespace("corpcor", quietly = TRUE)) {
      stop("Package 'corpcor' needed for nlshrink method")
    }
    Omega1 <- corpcor::invcov.shrink(G1, verbose = FALSE)
    Omega0 <- corpcor::invcov.shrink(G0, verbose = FALSE)
  } else if (method == "corpcor") {
    if (!requireNamespace("corpcor", quietly = TRUE)) {
      stop("Package 'corpcor' needed for corpcor method")
    }
    Omega1 <- corpcor::pseudoinverse(G1)
    Omega0 <- corpcor::pseudoinverse(G0)
  } else {
    stop("Method must be 'direct', 'nlshrink', or 'corpcor'")
  }

  return(list(Omega1 = Omega1, Omega0 = Omega0))
}


#' Compute basis projection matrices
#'
#' @param Z Basis matrix (n x k)
#' @param A treament weight (n x 1)
#' @param Omega1 Inverse Gram matrix for treatment group
#' @param Omega0 Inverse Gram matrix for control group
#'
#' @return List with B1 and B0 (projection matrices)
#' @export
compute_basis_matrix <- function(Z, A, Omega1, Omega0) {
  Z_weighted1 <- Z * A
  Z_weighted0<-  Z * (1-A)
  B1 <- Z %*% Omega1 %*% t(Z_weighted1)
  B0 <- Z %*% Omega0 %*% t(Z_weighted0)

  return(list(B1 = B1, B0 = B0))
}



#' Build E_j tensor structure
#'
#' Constructs a nested list representing the tensor structure [1, [1,2], ..., [j-1,j], j]
#' used in Higher-Order Influence Function (HOIF) calculations.
#'
#' @param j Integer greater than or equal to 2 specifying the dimension
#' @return A nested list with j+1 elements representing the tensor structure:
#'   - First element: scalar 1
#'   - Middle elements: vectors [1,2], [2,3], ..., [j-1,j]
#'   - Last element: scalar j
#' @export
build_Ej <- function(j) {
  if (j < 2) stop("j must be >= 2")

  E_j <- vector("list", j + 1)

  # First single index tensor: [1]
  E_j[[1]] <- 1

  # Middle B tensors: [1,2], [2,3], ..., [j-1,j]
  for (k in 1:(j-1)) {
    E_j[[k + 1]] <- c(k, k + 1)
  }

  # Last single index tensor: [j]
  E_j[[j + 1]] <- j

  E_j
}

#' Compute HOIF Estimators for ATE
#'
#' @param residuals A list containing the computed residuals: `R1`, `r1`, `R0`, and `r0`.
#' @param B_matrices A list containing the projection-like basis matrices: `B1` and `B0`.
#' @param m Integer. The maximum order of the HOIF estimator.
#' @param backend Character. The computation backend for the `ustat` package;
#'        either "torch" (default) or "numpy".
#' @param pure_R_code Logical. Whether to use a native R implementation.
#'        This serves as a fallback when the `reticulate` environment (Python)
#'        encounters configuration issues. Note: The pure R implementation
#'        only supports up to the 6th order ($m = 6$).
#'
#' @return A list of HOIF estimators (ATE, HOIF, IIFF) for orders $l = 2, \dots, m$.
#'
#' @importFrom ustats ustat
#' @export
compute_hoif_estimators <- function(residuals, B_matrices, m = 7, backend = "torch", pure_R_code = FALSE) {

  # Initialize storage
  U1 <- numeric(m)
  U0 <- numeric(m)
  IIFF1 <- numeric(m)
  IIFF0 <- numeric(m)
  HOIF1 <- numeric(m)
  HOIF0 <- numeric(m)
  ATE <- numeric(m)

  # Extract components
  R1 <- residuals$R1
  r1 <- residuals$r1
  R0 <- residuals$R0
  r0 <- residuals$r0
  B1 <- B_matrices$B1
  B0 <- B_matrices$B0

  if (!pure_R_code) {
    # Compute U-statistics for each order j = 2 to m using Python backend
    for (j in 2:m) {
      # Construct tensor list T_j^a
      # T_j^a = list(R^a, B^a, B^a, ..., B^a (j-1 times), r^a)
      T_j_1 <- c(list(r1), rep(list(B1), j - 1), list(R1))
      T_j_0 <- c(list(r0), rep(list(B0), j - 1), list(R0))

      # Construct expression E_j^a directly as nested list
      # For j tensors: [1, [1,2], [2,3], ..., [j-1,j], j]
      # This represents: R_{i_1} * B_{i_1,i_2} * B_{i_2,i_3} * ... * B_{i_{j-1},i_{j}} * r_{i_j}
      E_j <- build_Ej(j)

      # Compute U-statistics using ustat function
      U1[j] <- (-1)^j * ustats::ustat(tensors = T_j_1, expression = E_j,
                              backend = backend, average = TRUE)
      U0[j] <- (-1)^j * ustats::ustat(tensors = T_j_0, expression = E_j,
                              backend = backend, average = TRUE)
    }
  } else {
    # Compute U-statistics using pure R fallback (up to 6th order)
    U_list_1 <- calculate_u_statistics_pure_r_six(Vector_1 = r1, Vector_2 = R1,
                                           A1 = B1, A2 = B1, A3 = B1, A4 = B1, A5 = B1)
    U_list_0 <- calculate_u_statistics_pure_r_six(Vector_1 = r0, Vector_2 = R0,
                                           A1 = B0, A2 = B0, A3 = B0, A4 = B0, A5 = B0)
    for (j in 2:m) {
      U1[j] <- (-1)^j * U_list_1[[j - 1]]
      U0[j] <- (-1)^j * U_list_0[[j - 1]]
    }
  }
  # Compute IIFF and HOIF for each order l = 2 to m
  for (l in 2:m) {
    # IIFF_l = sum_{j=2}^l C_j^l * U_j
    # where C_j^l = choose(l-2, l-j)
    for (j in 2:l) {
      C_jl <- choose(l - 2, l - j)
      IIFF1[l] <- IIFF1[l] + C_jl * U1[j]
      IIFF0[l] <- IIFF0[l] + C_jl * U0[j]
    }

    # HOIF_l = sum_{j=2}^l IIFF_j
    HOIF1[l] <- sum(IIFF1[2:l])
    HOIF0[l] <- sum(IIFF0[2:l])

    # ATE_l = HOIF_l^1 - HOIF_l^0
    ATE[l] <- HOIF1[l] - HOIF0[l]
  }

  # Return results for orders 2 to m
  return(list(
    ATE = ATE[2:m],
    HOIF1 = HOIF1[2:m],
    HOIF0 = HOIF0[2:m],
    IIFF1 = IIFF1[2:m],
    IIFF0 = IIFF0[2:m],
    orders = 2:m
  ))
}



#' Main function: HOIF estimators for ATE with optional sample splitting
#'
#' @param X Covariate matrix (n x p)
#' @param A Treatment vector (n x 1)
#' @param Y Outcome vector (n x 1)
#' @param mu1 Predicted outcomes under treatment
#' @param mu0 Predicted outcomes under control
#' @param pi Predicted propensity scores
#' @param transform_method Character: method to transform covariates before
#'   constructing basis functions.
#'   - "splines": use basis splines expansion
#'   - "fourier": use Fourier basis expansion
#'   - "none": no transformation (use raw covariates)
#'
#' @param basis_dim Integer: number of basis functions to generate when using
#'   "splines" or "fourier" transformations. Higher values provide more flexible
#'   approximations but may increase variance.
#'
#' @param inverse_method Character: regularization method for Gram matrix inversion.
#'   - "direct": direct Moore-Penrose pseudoinverse (no regularization)
#'   - "nlshrink": nonlinear shrinkage estimator (Ledoit-Wolf type)
#'   - "corpcor": shrinkage via the corpcor package (for high-dimensional settings)
#' @param m Maximum order for HOIF
#' @param sample_split Logical: whether to use sample splitting.
#'   If `TRUE`, the data is split: one part for estimating the inverse
#'   Gram matrix, and the other for estimation. If `FALSE`, it corresponds
#'   to the sHOIF case (without sample splitting).
#' @param n_folds Number of folds for sample splitting (if used)
#' @param backend Character: "torch" (default) or "numpy"
#' @param seed Random seed for reproducibility (for sample splitting)
#' @param ... Additional arguments passed to transform_covariates
#'
#' @return List with ATE, HOIF, and IIFF estimates
#' @seealso \code{\link{compute_HOIF_test}}, which provides a brute-force
#'   implementation used internally for correctness checks on small datasets.
#' @export
hoif_ate <- function(X, A, Y, mu1, mu0, pi,
                     transform_method = "splines",
                     basis_dim,
                     inverse_method = "direct",
                     m = 7,
                     sample_split = FALSE,
                     n_folds = 2,
                     backend = "torch",
                     seed = NULL,
                     pure_R_code = FALSE,
                     ...) {

  n <- length(Y)

  # Check for order limit in pure R mode
  if (pure_R_code) {
    if (missing(m)) {
      m <- 6
      warning("Using default order m = 6 for pure R mode. The maximum order supported in pure R is 6.")
    } else if (m > 6) {
      stop("When 'pure_R_code' is TRUE, the maximum order supported is 6. If you want order > 6, please set 'pure_R_code' be FALSE then use powerful python. ")
    } else if (m < 6) {
      warning("You requested order ", m, " with 'pure_R_code' = TRUE. ",
              "The pure R implementation computes all orders up to 6 at once. ",
              "Therefore, the calculation time is the same as for order 6, ",
              "and results up to 6th order will be returned.")
      m <- 6
    }
  }

  # Step 1: Transform covariates (done on full data)
  Z <- transform_covariates(X, method = transform_method, basis_dim = basis_dim, ...)

  # Step 2: Compute residuals (done on full data)
  residuals <- compute_residuals(A, Y, mu1, mu0, pi)

  if (!sample_split) {
    # No sample splitting: standard procedure

    # Step 3: Compute inverse Gram matrices
    Omega <- compute_gram_inverse(Z, A, method = inverse_method)

    # Step 4: Compute basis matrices
    B_matrices <- compute_basis_matrix(Z, A, Omega$Omega1, Omega$Omega0)

    # Step 5: Compute HOIF estimators
    results <- compute_hoif_estimators(residuals, B_matrices, m, backend, pure_R_code)

  } else {
    # Sample splitting (cross-fitting)

    # Set seed for fold creation if provided
    if (!is.null(seed)) {
      set.seed(seed)
    }

    # Create fold indices
    fold_indices <- sample(rep(1:n_folds, length.out = n))

    # Storage for fold-specific estimates
    ATE_folds <- matrix(0, nrow = n_folds, ncol = m - 1)
    HOIF1_folds <- matrix(0, nrow = n_folds, ncol = m - 1)
    HOIF0_folds <- matrix(0, nrow = n_folds, ncol = m - 1)
    IIFF1_folds <- matrix(0, nrow = n_folds, ncol = m - 1)
    IIFF0_folds <- matrix(0, nrow = n_folds, ncol = m - 1)

    for (j in 1:n_folds) {
      # Define fold indices
      I_j <- which(fold_indices == j)
      I_not_j <- which(fold_indices != j)

      # Step 3: Compute Omega on training set (not I_j)
      Z_train <- Z[I_not_j, , drop = FALSE]
      A_train <- A[I_not_j]
      Omega_j <- compute_gram_inverse(Z_train, A_train, method = inverse_method)

      # Step 4: Compute basis matrices on test set (I_j)
      Z_test <- Z[I_j, , drop = FALSE]
      A_test <- A[I_j]
      B_matrices_j <- compute_basis_matrix(Z_test, A_test, Omega_j$Omega1, Omega_j$Omega0)

      # Step 5: Compute HOIF on test set
      residuals_j <- list(
        R1 = residuals$R1[I_j],
        r1 = residuals$r1[I_j],
        R0 = residuals$R0[I_j],
        r0 = residuals$r0[I_j]
      )

      results_j <- compute_hoif_estimators(residuals_j, B_matrices_j, m, backend, pure_R_code)

      # Store results
      ATE_folds[j, ] <- results_j$ATE
      HOIF1_folds[j, ] <- results_j$HOIF1
      HOIF0_folds[j, ] <- results_j$HOIF0
      IIFF1_folds[j, ] <- results_j$IIFF1
      IIFF0_folds[j, ] <- results_j$IIFF0
    }

    # Average across folds
    results <- list(
      ATE = colMeans(ATE_folds),
      HOIF1 = colMeans(HOIF1_folds),
      HOIF0 = colMeans(HOIF0_folds),
      IIFF1 = colMeans(IIFF1_folds),
      IIFF0 = colMeans(IIFF0_folds),
      orders = 2:m
    )
}

  # Add convergence plot data
  results$convergence_data <- data.frame(
    order = results$orders,
    ATE = results$ATE
  )

  class(results) <- "hoif_ate"
  return(results)
}


#' Print method for hoif_ate objects
#'
#' @param x Object of class hoif_ate
#' @param ... Additional arguments (unused)
#'
#' @export
print.hoif_ate <- function(x, ...) {
  cat("HOIF Estimators for Average Treatment Effect\n")
  cat("=============================================\n\n")

  cat("Estimates by order:\n")
  print(data.frame(
    Order = x$orders,
    ATE = round(x$ATE, 4),
    HOIF1 = round(x$HOIF1, 4),
    HOIF0 = round(x$HOIF0, 4)
  ))

  cat("\nFinal ATE estimate (highest order):", round(tail(x$ATE, 1), 4), "\n")
}


#' Plot convergence of ATE estimates
#'
#' @param x Object of class hoif_ate
#' @param ... Additional arguments passed to plot
#'
#' @export
plot.hoif_ate <- function(x, ...) {
  plot(x$orders, x$ATE, type = "b", pch = 19,
       xlab = "Order", ylab = "ATE Estimate",
       main = "Convergence of HOIF-ATE Estimator",
       ...)
  abline(h = tail(x$ATE, 1), lty = 2, col = "red")
  grid()
}


#--------------------------U-stat pure in R--------------------
#' Compute U-statistics HOIF-type from Order 2 to 6 in pure R code.
#'
#' This function serves as a core computational component in higher-order
#' influence function (HOIF) estimators in pure R code.
#'
#' Internally, the function constructs kernel matrices for orders 2 through 6
#' using recursive matrix operations and removes diagonal contributions to
#' ensure degenerate U-statistics.
#'
#' @param Vector_1 Numeric column vector of length \eqn{n}.
#' @param Vector_2 Numeric column vector of length \eqn{n}.
#' @param A1,A2,A3,A4,A5 Numeric \eqn{n \times n} kernel matrices of the same dimension.
#'
#' @return A named list containing numeric U-statistic estimates:
#' \describe{
#'   \item{U2}{Second-order U-statistic}
#'   \item{U3}{Third-order U-statistic}
#'   \item{U4}{Fourth-order U-statistic}
#'   \item{U5}{Fifth-order U-statistic}
#'   \item{U6}{Sixth-order U-statistic}
#' }
#'
#'@importFrom SMUT eigenMapMatMult
#' @details
#' All diagonal elements of intermediate kernel matrices are removed to avoid
#' self-interactions. Matrix multiplications are performed via
#' `eigenMapMatMult()` and element-wise products via `hadamard()`.
#' The exact formula of the output is:
#' \deqn{
#' \mathbb{U}_{n,m} = \frac{1}{\binom{n}{m} m!}
#' \sum_{i_1 \ne \cdots \ne i_m} Vector_1[i_1] \cdot A1[i_1,i_2] \cdot A1[i_2,i_3] \cdots A1[i_{m-1},i_{m}] \cdot Vector_2[i_{m}]
#' }
#'
#' @export
calculate_u_statistics_pure_r_six <- function(Vector_1, Vector_2, A1, A2 , A3, A4, A5) {
  # Ensure input matrices are square and of the same size
  n <- nrow(A1)
  if (!(nrow(A2) == n &&
        nrow(A3) == n && nrow(A4) == n && nrow(A5) == n)) {
    stop("All input matrices must be square and of the same size")
  }

  # Precompute no_diag matrices to avoid recomputation
  no_diag_A1 <- no_diag(A1)
  no_diag_A2 <- no_diag(A2)
  no_diag_A3 <- no_diag(A3)
  no_diag_A4 <- no_diag(A4)
  no_diag_A5 <- no_diag(A5)

  Ker_1234 <-  calculate_u_Ker_5(A1, A2, A3, A4)
  Ker5 <- Ker_1234$Ker5
  Ker4 <- Ker_1234$Ker4
  Ker3 <- Ker_1234$Ker3
  Ker2 <- Ker_1234$Ker2

  # Calculate Sub matrices for 6th order
  # Sub matrix for i2 = i6
  Ker_2345 <- calculate_u_Ker_5(A2, A3, A4, A5)
  Ker5_2345 <- Ker_2345$Ker5
  Ker4_234 <- Ker_2345$Ker4

  Ker_345 <- calculate_u_Ker_4(A3, A4, A5)
  Ker4_345 <- Ker_345$Ker4

  Sub_6_2 <- eigenMapMatMult(no_diag_A1, diag_Mat(Ker5_2345)) -
    hadamard(hadamard(no_diag_A1, t(no_diag_A2)), Ker4_345) -
    hadamard(hadamard(no_diag_A1, t(no_diag(Ker4_234))), no_diag_A5) -
    hadamard(hadamard(no_diag_A1, eigenMapMatMult(t(no_diag_A3), t(no_diag_A2))),eigenMapMatMult(no_diag_A4, no_diag_A5)) +
    hadamard(no_diag_A1, eigenMapMatMult(hadamard(no_diag_A4, t(no_diag_A3)), hadamard(no_diag_A5, t(no_diag_A2))))

  # Sub matrix for i3 = i6
  Sub_6_3 <-     eigenMapMatMult(eigenMapMatMult(no_diag_A1, no_diag_A2), diag_Mat(Ker4_345)) -
    hadamard(hadamard(eigenMapMatMult(no_diag_A1, no_diag_A2), t(no_diag_A3)),eigenMapMatMult(no_diag_A4, no_diag_A5)) -
    hadamard(  hadamard(eigenMapMatMult(no_diag_A1, no_diag_A2), t(eigenMapMatMult(  no_diag_A3, no_diag_A4  ))), no_diag_A5)-
    eigenMapMatMult(no_diag_A1, hadamard(  hadamard(no_diag_A2, t(no_diag_A3)),eigenMapMatMult(no_diag_A4, no_diag_A5))) +
    hadamard(eigenMapMatMult(hadamard(no_diag_A1,t(no_diag_A4)), hadamard(no_diag_A2, t(no_diag_A3))), no_diag_A5) -
    eigenMapMatMult(no_diag_A1, hadamard(    hadamard(no_diag_A2, no_diag_A5), t(eigenMapMatMult(no_diag_A3, no_diag_A4))  )) +
    hadamard(t(no_diag_A3), eigenMapMatMult(hadamard(no_diag_A1, no_diag_A4), hadamard(no_diag_A2, no_diag_A5) ))


  # Sub matrix for i4 = i6
  Sub_6_4 <-   eigenMapMatMult(Ker4, diag_Mat(eigenMapMatMult(no_diag_A4, no_diag_A5)))-
    hadamard(Ker4, hadamard(t(no_diag_A4), no_diag_A5)) -
    eigenMapMatMult(no_diag_A1, hadamard( eigenMapMatMult(no_diag_A2, no_diag_A3), hadamard(t(no_diag_A4), no_diag_A5)) ) +
    hadamard(no_diag_A3, eigenMapMatMult(hadamard(no_diag_A1, t(no_diag_A2)), hadamard(t(no_diag_A4), no_diag_A5))) -
    eigenMapMatMult(no_diag_A1,no_diag(eigenMapMatMult(no_diag_A2, hadamard(hadamard( no_diag_A3, t(no_diag_A4)), no_diag_A5)))) +
    eigenMapMatMult( diag_col_sum(hadamard(t(no_diag_A1),no_diag_A2)), hadamard(hadamard(no_diag_A3, t(no_diag_A4)),no_diag_A5)) -
    hadamard(hadamard(hadamard(hadamard(no_diag_A1, t(no_diag_A2)), no_diag_A3), t(no_diag_A4)), no_diag_A5)

  # Kernel matrix for 6th order
  Ker6 <- eigenMapMatMult(no_diag(Ker5), no_diag_A5) - Sub_6_2 - Sub_6_3 - Sub_6_4

  # Calculate U-statistics
  # For each order, we use the formula from the document
  Ker2 <- no_diag(Ker2)
  Ker3 <- no_diag(Ker3)
  Ker4 <- no_diag(Ker4)
  Ker5 <- no_diag(Ker5)
  Ker6 <- no_diag(Ker6)

  U2_all <- eigenMapMatMult(eigenMapMatMult(t(Vector_1), Ker2), Vector_2) / (n * (n-1))
  U3_all <- eigenMapMatMult(eigenMapMatMult(t(Vector_1), Ker3), Vector_2) / (n * (n-1) * (n-2))
  U4_all <- eigenMapMatMult(eigenMapMatMult(t(Vector_1), Ker4), Vector_2) / (n * (n-1) * (n-2) * (n-3))
  U5_all <- eigenMapMatMult(eigenMapMatMult(t(Vector_1), Ker5), Vector_2) / (n * (n-1) * (n-2) * (n-3) * (n-4))
  U6_all <- eigenMapMatMult(eigenMapMatMult(t(Vector_1), Ker6), Vector_2) / (n * (n-1) * (n-2) * (n-3) * (n-4) * (n-5) )

  results_list <- list(
    U2 = as.numeric(U2_all),
    U3 = as.numeric(U3_all),
    U4 = as.numeric(U4_all),
    U5 = as.numeric(U5_all),
    U6 = as.numeric(U6_all)
  )
  # Return a list of U-statistics
  return( results_list)
}

#' Remove Diagonal Elements of a Matrix
#'
#' Sets the diagonal entries of a matrix to zero. This is used to eliminate
#' self-interactions when constructing U-statistic kernel matrices.
#'
#' @param mat A numeric square matrix.
#'
#' @return The same matrix with its diagonal set to zero.
#' @keywords internal
no_diag <- function(mat) {
  diag(mat) <- 0
  return(mat)
}

#' Hadamard (Element-wise) Matrix Product
#'
#' Computes the element-wise (Hadamard) product of two matrices of the same size.
#'
#' @param mat1,mat2 Numeric matrices of identical dimensions.
#'
#' @return A matrix containing the element-wise product.
#' @keywords internal
hadamard <- function(mat1, mat2) {
  return(mat1 * mat2)
}

#' Diagonal Matrix of Column Sums
#'
#' Creates a diagonal matrix whose diagonal elements are the column sums
#' of the input matrix.
#'
#' @param mat A numeric matrix.
#'
#' @return A diagonal matrix.
#' @keywords internal
diag_col_sum <- function(mat) {
  diag(colSums(mat))
}
#' Extract Diagonal as a Diagonal Matrix
#'
#' Returns a diagonal matrix containing only the diagonal elements of
#' the input matrix.
#'
#' @param mat A numeric square matrix.
#'
#' @return A diagonal matrix formed from \code{diag(mat)}.
#' @keywords internal
diag_Mat <- function(mat) {
  diag(diag(mat))
}
#' Construct Kernel Matrices up to 4th Order
#'
#' Builds second-, third-, and fourth-order kernel matrices used in the
#' recursive construction of higher-order U-statistics.
#'
#' @param A1,A2,A3 Numeric square matrices of the same dimension.
#' @importFrom SMUT eigenMapMatMult
#' @return A list containing kernel matrices \code{Ker2}, \code{Ker3}, and \code{Ker4}.
#' @keywords internal
calculate_u_Ker_4 <- function(A1, A2, A3) {
  # Ensure input matrices are square and of the same size
  n <- nrow(A1)
  if (!(nrow(A2) == n &&
        nrow(A3) == n)) {
    stop("All input matrices must be square and of the same size")
  }

  # Precompute no_diag matrices to avoid recomputation
  no_diag_A1 <- no_diag(A1)
  no_diag_A2 <- no_diag(A2)
  no_diag_A3 <- no_diag(A3)


  # Kernel matrix for 2nd order
  Ker2 <- A1

  # Kernel matrix for 3rd order
  Ker3 <- eigenMapMatMult(no_diag_A1, no_diag_A2)

  # Kernel matrix for 4th order
  Ker4 <- eigenMapMatMult(no_diag(Ker3), no_diag_A3) -
    eigenMapMatMult(A1, diag_col_sum(t(no_diag_A2) * no_diag_A3)) +
    hadamard(hadamard(A1, t(A2)), A3)
  diag(Ker4) <- diag(eigenMapMatMult(Ker3, no_diag_A3))
  return(list(
    Ker2 = Ker2,
    Ker3 = Ker3,
    Ker4 = Ker4
  ))
}
#' Construct Kernel Matrices up to 5th Order
#'
#' Extends the recursive kernel construction to fifth order based on four
#' input kernel matrices. This function builds upon
#' \code{calculate_u_Ker_4()}.
#'
#' @param A1,A2,A3,A4 Numeric square matrices of the same dimension.
#' @importFrom SMUT eigenMapMatMult
#' @return A list containing kernel matrices \code{Ker2}, \code{Ker3},
#' \code{Ker4}, and \code{Ker5}.
#' @keywords internal
calculate_u_Ker_5 <- function(A1, A2, A3, A4) {
  # Ensure input matrices are square and of the same size
  n <- nrow(A1)
  if (!(nrow(A2) == n &&
        nrow(A3) == n && nrow(A4) == n)) {
    stop("All input matrices must be square and of the same size")
  }

  # Precompute no_diag matrices to avoid recomputation
  no_diag_A1 <- no_diag(A1)
  no_diag_A2 <- no_diag(A2)
  no_diag_A3 <- no_diag(A3)
  no_diag_A4 <- no_diag(A4)

  Ker_123 <- calculate_u_Ker_4(A1,A2,A3)
  Ker4 <-   Ker_123$Ker4
  Ker3 <-   Ker_123$Ker3
  Ker2 <-   Ker_123$Ker2

  # Kernel matrix for 5th order (more complex)
  Ker5 <- eigenMapMatMult(no_diag(Ker4), no_diag_A4) -
    eigenMapMatMult(no_diag_A1, diag(diag(eigenMapMatMult(eigenMapMatMult(no_diag_A2, no_diag_A3), no_diag_A4)))) +
    hadamard(hadamard(no_diag_A1, t(no_diag_A2)), eigenMapMatMult(no_diag_A3, no_diag_A4)) +
    hadamard(hadamard(no_diag_A1, no_diag_A4), eigenMapMatMult(t(no_diag_A3), t(no_diag_A2))) -
    eigenMapMatMult(no_diag(eigenMapMatMult(no_diag_A1, no_diag_A2)), diag(diag(eigenMapMatMult(no_diag_A3, no_diag_A4)))) +
    eigenMapMatMult(no_diag_A1, no_diag(hadamard(hadamard(no_diag_A2, t(no_diag_A3)), no_diag_A4))) +
    hadamard(eigenMapMatMult(no_diag_A1, no_diag_A2), hadamard(t(no_diag_A3), no_diag_A4))
  diag(Ker5) <- diag(eigenMapMatMult(Ker4, no_diag_A4))
  return(list(
    Ker2 = Ker2,
    Ker3 = Ker3,
    Ker4 = Ker4,
    Ker5 = Ker5
  ))
}


#----------------------------test function brute force calculation (only for double check the computation test) ---------------------
#' Generate All Ordered m-Tuples of Distinct Indices (Brute Force)
#'
#' Internal helper function used for brute-force validation of higher-order
#' influence function (HOIF) calculations. It generates all ordered
#' \eqn{m}-permutations of a given index set.
#'
#' ⚠️ **Warning:** This function has factorial computational complexity and
#' should only be used for very small sample sizes as part of debugging or
#' verification routines.
#'
#' @param indices Integer vector of indices.
#' @param m Integer order of the permutation.
#'
#' @return A matrix where each row is an ordered \eqn{m}-tuple of distinct indices.
#'
#' @keywords internal
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

#' Generate All Permutations of a Vector (Recursive, Brute Force)
#'
#' Internal recursive helper used by \code{generate_permutations()} to enumerate
#' all permutations of a vector.
#'
#' ⚠️ Extremely computationally expensive for vectors longer than ~8 elements.
#' Only intended for internal brute-force validation code.
#'
#' @param x A vector.
#'
#' @return A list of vectors, each being a permutation of \code{x}.
#'
#' @keywords internal
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
#' Brute-Force Sequence of Higher-Order Influence Function Estimates (Test Only)
#'
#' Computes a sequence of higher-order influence function (HOIF) estimates
#' from order 2 up to order \code{m} using a brute-force implementation.
#' This function is **only intended for double-checking correctness** of the
#' main fast implementation.
#'
#' ⚠️ The computation scales combinatorially with both sample size and order
#' and becomes infeasible beyond very small datasets.
#'
#' @param X Covariate matrix (n × p).
#' @param A Binary treatment vector of length n.
#' @param Y Outcome vector of length n.
#' @param mu1 Estimated outcome regression under treatment.
#' @param mu0 Estimated outcome regression under control.
#' @param pi Estimated propensity scores.
#' @param m Maximum order of HOIF to compute.
#' @param sample_splitting Logical or integer; whether to use sample splitting.
#' @param n_folds Number of folds for sample splitting.
#' @param seed Random seed for fold assignment.
#'
#' @return A list containing cumulative HOIF and incremental IIFF terms for
#' both treatment arms.
#'
#' @details
#' This routine repeatedly calls \code{compute_HOIF_test()} and accumulates
#' incremental influence function terms. It is not optimized and should never
#' be used in production or large-sample simulations.
#'
#' @keywords internal
compute_HOIF_sequence_test <- function(X, A, Y, mu1, mu0, pi, m, sample_splitting = 0, n_folds = 2, seed) {
  HOIF_1_m <- numeric(m-1)
  HOIF_0_m <- numeric(m-1)
  IIFF_1_m <- numeric(m-1)
  IIFF_0_m <- numeric(m-1)
  for (i in 2:m)
  {
    results <- compute_HOIF_test(X, A, Y, mu1, mu0, pi, i, sample_splitting, n_folds,seed )
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
#' Brute-Force Higher-Order Influence Function Estimator (Test Version)
#'
#' Computes the order-\code{m} higher-order influence function (HOIF) term
#' using an explicit enumeration of all ordered index tuples.
#'
#' 🚨 **This implementation is intentionally naive and is provided solely
#' for validation and debugging purposes.** It should only be used on very
#' small datasets to verify the correctness of optimized implementations.
#'
#' The computational complexity grows factorially with both sample size and
#' order \code{m}.
#'
#' @param X Covariate matrix (n × p).
#' @param A Binary treatment indicator vector.
#' @param Y Outcome vector.
#' @param mu1 Estimated outcome regression under treatment.
#' @param mu0 Estimated outcome regression under control.
#' @param pi Estimated propensity scores.
#' @param m Order of the HOIF term.
#' @param sample_splitting Whether to use K-fold sample splitting (0 = no).
#' @param n_folds Number of folds for sample splitting.
#' @param seed Random seed for reproducibility.
#'
#' @return A list with elements:
#' \describe{
#'   \item{IIFF_1_m}{Order-\code{m} influence function term for treated units}
#'   \item{IIFF_0_m}{Order-\code{m} influence function term for control units}
#' }
#'
#' @details
#' This function directly implements the combinatorial definition of the HOIF
#' using full enumeration of index permutations. Matrix inverses are computed
#' explicitly and no numerical stabilization is included.
#'
#' This function is part of the package's **internal test infrastructure** and
#' should not be exported or used in applied analysis.
#'
#' @keywords internal
compute_HOIF_test <- function(X, A, Y, mu1, mu0, pi, m, sample_splitting = 0, n_folds = 2, seed) {
  # X: n x p covariate matrix
  # A: n-vector of binary treatment (0 or 1)
  # Y: n-vector of outcomes
  # mu1: n-vector of estimated mu(1, X_i)
  # mu0: n-vector of estimated mu(0, X_i)
  # pi: n-vector of estimated propensity scores
  # m: order of the HOIF statistic
  # sample_splitting: 0 for no splitting, 1 for K-fold splitting
  # n_folds: number of folds (only used when sample_splitting = 1)

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
      term <- r_a_est[idx[1]] * (t(Z_est[idx[1], , drop = FALSE]) %*% Omega_a)

      # Middle products for s = 2 to m-1
      if (m >= 3) {
        for (s in 2:(m-1)) {
          Q_a_is <- s_a_est[idx[s]] * Z_est[idx[s], , drop = FALSE] %*% t(Z_est[idx[s], , drop = FALSE])



          term <- term %*% (Q_a_is - Sigma_a) %*% Omega_a
          # term <- term %*% (Q_a_is ) %*% Omega_a
        }
      }

      # Final term: Z_{i_m} * r^a_{i_m}
      term <- term %*% Z_est[idx[m], , drop = FALSE] * R_a_est[idx[m]] *  s_a_est[idx[m]]

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

    IIFF_1_m <- numeric(n_folds)
    IIFF_0_m <- numeric(n_folds)

    for (j in 1:n_folds) {
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
