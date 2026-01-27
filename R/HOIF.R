#' HOIF Estimators for Average Treatment Effect
#'
#' This file implements the Higher Order Influence Function (HOIF) estimators
#' for Average Treatment Effect (ATE) estimation with nuisance functions.
#'
#' @author Xingyu Chen




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
  R1 <- A * (Y - mu1)
  r1 <- 1 - A / pi

  # For a = 0
  R0 <- (1 - A) * (Y - mu0)
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
#' @param Omega1 Inverse Gram matrix for treatment group
#' @param Omega0 Inverse Gram matrix for control group
#'
#' @return List with B1 and B0 (projection matrices)
#' @export
compute_basis_matrix <- function(Z, Omega1, Omega0) {
  B1 <- Z %*% Omega1 %*% t(Z)
  B0 <- Z %*% Omega0 %*% t(Z)

  return(list(B1 = B1, B0 = B0))
}


#' Convert nested list expression to Einstein notation
#'
#' Converts a nested list expression like [[1,2],[2,3],[3,4]] to Einstein
#' notation like "ab,bc,cd->".
#'
#' @param expr_list A list of integer vectors, each of length 2
#'
#' @return Character string in Einstein notation
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' expr_list_to_einstein([[1,2],[2,3],[3,4]])  # Returns "ab,bc,cd->"
#' }
expr_list_to_einstein <- function(expr_list) {
  # Step 1: Collect all unique indices (as characters for safe naming)
  all_indices <- sort(unique(unlist(expr_list)))

  if (length(all_indices) > 26) {
    stop("Too many unique indices (>26); Einstein notation limited to a-z")
  }

  # Step 2: Map each index to a letter: 1->a, 2->b, ..., 26->z
  # Use as.character() to avoid numeric name issues
  index_to_letter <- setNames(letters[seq_along(all_indices)],
                              as.character(all_indices))

  # Step 3: Convert each group (e.g., c(1,2,3) -> "abc")
  terms <- sapply(expr_list, function(indices) {
    # Ensure indices is a vector (even if length 1)
    idx_chars <- as.character(indices)
    letters_vec <- index_to_letter[idx_chars]

    # Safety check: any missing?
    if (any(is.na(letters_vec))) {
      stop("Found unmapped index in: ", paste(indices, collapse = ", "))
    }

    paste0(letters_vec, collapse = "")
  })

  # Step 4: Join with commas and append "->"
  paste0(paste(terms, collapse = ","), "->")
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
      T_j_1 <- c(list(R1), rep(list(B1), j - 1), list(r1))
      T_j_0 <- c(list(R0), rep(list(B0), j - 1), list(r0))

      # Construct expression E_j^a directly as nested list
      # For j tensors: [1, [1,2], [2,3], ..., [j-1,j], j]
      # This represents: R_i * B_{i,i1} * B_{i1,i2} * ... * r_{ij}
      E_j <- build_Ej(j)

      # Compute U-statistics using ustat function
      U1[j] <- (-1)^j * ustats::ustat(tensors = T_j_1, expression = E_j,
                              backend = backend, average = TRUE)
      U0[j] <- (-1)^j * ustats::ustat(tensors = T_j_0, expression = E_j,
                              backend = backend, average = TRUE)
    }
  } else {
    # Compute U-statistics using pure R fallback (up to 6th order)
    U_list_1 <- calculate_u_statistics_six(Vector_1 = R1, Vector_2 = r1,
                                           A1 = B1, A2 = B1, A3 = B1, A4 = B1, A5 = B1)
    U_list_0 <- calculate_u_statistics_six(Vector_1 = R0, Vector_2 = r0,
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
    B_matrices <- compute_basis_matrix(Z, Omega$Omega1, Omega$Omega0)

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
      B_matrices_j <- compute_basis_matrix(Z_test, Omega_j$Omega1, Omega_j$Omega0)

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
calculate_u_statistics_six <- function(Vector_1, Vector_2, A1, A2 , A3, A4, A5) {
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

