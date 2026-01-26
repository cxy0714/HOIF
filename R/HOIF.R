#' HOIF Estimators for Average Treatment Effect
#'
#' This file implements the Higher Order Influence Function (HOIF) estimators
#' for Average Treatment Effect (ATE) estimation with nuisance functions.
#'
#' @author Xingyu Chen
#' @date 2026-01-23

# Required packages
# install.packages(c("splines", "corpcor"))

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
#' expr_list_to_einstein(list(c(1,2), c(2,3), c(3,4)))  # Returns "ab,bc,cd->"
#' }
expr_list_to_einstein <- function(expr_list) {
  # Collect and SORT indices
  all_indices <- sort(unique(as.integer(unlist(expr_list))))
  n_unique <- length(all_indices)

  if (n_unique > 26) {
    stop("Too many unique indices (>26)")
  }

  index_to_letter <- setNames(letters[seq_len(n_unique)], all_indices)

  # Convert each tensor
  terms <- vapply(expr_list, function(idx) {
    idx <- as.integer(idx)

    if (length(idx) == 1) {
      index_to_letter[as.character(idx)]

    } else if (length(idx) == 2) {
      paste0(index_to_letter[as.character(idx[1])],
             index_to_letter[as.character(idx[2])])

    } else {
      stop("Only rank-1 or rank-2 tensors supported")
    }
  }, character(1))

  paste0(paste(terms, collapse = ","), "->")
}



#' Compute a Higher-Order U-Statistic
#'
#' Computes a higher-order U-statistic given pre-computed kernel matrices
#' using either an Einstein summation expression or nested list notation.
#' This function calls Python's u-stats package via reticulate.
#'
#' **Note**: This function requires Python and the u-stats package to be installed.
#' Run \code{setup_hoif()} to install the required dependencies.
#'
#' @param tensors A list of numeric matrices/vectors of equal dimensions.
#' @param expression Either a character string (Einstein notation like "ab,bc->")
#'        or a nested list (like list(c(1,2), c(2,3))). Can also accept the special
#'        format list(1, list(1,2), ..., j) from compute_hoif_estimators.
#' @param backend Character string, either "numpy" or "torch" (default).
#' @param average Logical; whether to return the averaged U-statistic (default TRUE).
#'
#' @return A numeric scalar.
#'
#' @examples
#' \dontrun{
#' # First, setup the Python environment
#' setup_hoif()
#'
#' # Then use the function
#' v1 c- runif(100)
#' H1 <- matrix(runif(100), 10, 10)
#' H2 <- matrix(runif(100), 10, 10)
#' ustat(list(H1, H2), "ab,bc->")
#' # Or equivalently:
#' ustat(list(v1, H1, H2),"a,ab,bc->")
#' }
#'
#' @export
ustat <- function(tensors,
                  expression,
                  backend = c("torch", "numpy"),
                  average = TRUE) {

  if (!check_python_env()) {
    stop(
      "Python environment not properly configured.\n",
      "Please run: setup_hoif()\n",
      "Or check setup status with: check_hoif_setup()",
      call. = FALSE
    )
  }

  backend <- match.arg(backend)

  if (backend == "torch" && !reticulate::py_module_available("torch")) {
    warning(
      "Torch backend not available; falling back to numpy.\n",
      "Install torch with: reticulate::py_install('torch')",
      call. = FALSE
    )
    backend <- "numpy"
  }

  ustats <- tryCatch({
    reticulate::import("u_stats", delay_load = TRUE)
  }, error = function(e) {
    stop(
      "Failed to import u_stats module.\n",
      "Please run: setup_hoif()\n",
      "Error: ", e$message,
      call. = FALSE
    )
  })

  ustats$set_backend(backend)
  np <- reticulate::import("numpy", convert = FALSE)

  # ==========================================================
  # 🔁 Expression auto-conversion
  # ==========================================================
  if (is.list(expression)) {

    # 支持 build_Ej() 生成的格式: list(1, c(1,2), ..., j)
    valid_structure <- all(vapply(expression, function(x) {
      is.numeric(x) && length(x) %in% c(1, 2)
    }, logical(1)))

    if (!valid_structure) {
      stop("Expression list must contain numeric vectors of length 1 or 2.",
           call. = FALSE)
    }

    expression <- expr_list_to_einstein(expression)

  } else if (!is.character(expression)) {
    stop("Expression must be either a character string or a nested list",
         call. = FALSE)
  }

  # ==========================================================
  # 🔁 Tensor auto-conversion
  # ==========================================================
  tensors <- lapply(tensors, function(x) {

    # ✅ 已经是 Python numpy / torch 张量 → 直接放行
    if (inherits(x, "python.builtin.object")) {
      return(x)
    }

    # ❌ 必须是 numeric
    if (!is.numeric(x)) {
      stop("All tensors must be numeric (vector, matrix, or array).",
           call. = FALSE)
    }

    # 🔹 R 向量（无 dim） → 1D numpy array
    if (is.null(dim(x))) {
      return(np$array(as.numeric(x), dtype = "float32"))
    }

    # 🔹 任意维 R array → 保持维度转 numpy
    return(np$array(x, dtype = "float32"))
  })


  # ==========================================================
  # 🚀 Call Python ustat
  # ==========================================================
  result <- tryCatch({
    ustats$ustat(
      tensors = tensors,
      expression = expression,
      average = average
    )
  }, error = function(e) {
    stop(
      "Error computing U-statistic: ", e$message, "\n",
      "Expression used: ", expression,
      call. = FALSE
    )
  })

  as.numeric(result)
}


#' Transform covariates to basis functions
#'
#' @param X Matrix of covariates (n x p)
#' @param method Character: "splines" or "fourier"
#' @param k Integer: dimension of basis expansion
#' @param degree Integer: degree for B-splines (default 3)
#' @param period Numeric: period for Fourier basis (default 1)
#'
#' @return Matrix Z (n x k) of transformed covariates with intercept
#' @export
transform_covariates <- function(X, method = "splines", k = 10,
                                 degree = 3, period = 1) {
  if (is.vector(X)) X <- matrix(X, ncol = 1)
  n <- nrow(X)
  p <- ncol(X)

  # Scale X to [0, 1] for each dimension
  X_scaled <- apply(X, 2, function(x) {
    (x - min(x)) / (max(x) - min(x) + 1e-10)
  })

  # Initialize with intercept
  Z <- matrix(1, nrow = n, ncol = 1)

  # Generate basis for each dimension separately
  if (method == "splines") {
    for (j in 1:p) {
      # B-splines basis
      knots <- seq(0, 1, length.out = k - degree + 1)
      basis_j <- splines::bs(X_scaled[, j], knots = knots[-c(1, length(knots))],
                             degree = degree, intercept = FALSE)
      Z <- cbind(Z, basis_j)
    }
  } else if (method == "fourier") {
    for (j in 1:p) {
      # Fourier basis: cos and sin terms
      freq <- 1:floor(k / 2)
      for (f in freq) {
        Z <- cbind(Z, cos(2 * pi * f * X_scaled[, j] / period))
        Z <- cbind(Z, sin(2 * pi * f * X_scaled[, j] / period))
      }
    }
  } else {
    stop("Method must be 'splines' or 'fourier'")
  }

  return(Z)
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
  if (requireNamespace("sumt", quietly = TRUE)) {
    # Use faster matrix multiplication
    Z_weighted1 <- Z * sqrt(s1)
    Z_weighted0 <- Z * sqrt(s0)
    G1 <- sumt::eigenMapMatMult(t(Z_weighted1), Z_weighted1) / n
    G0 <- sumt::eigenMapMatMult(t(Z_weighted0), Z_weighted0) / n
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



#' @param j Integer >= 2
#' @return Nested list representing [1,[1,2],...,[j-1,j],j]
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

#' Compute HOIF estimators for ATE
#'
#' @param residuals List with R1, r1, R0, r0
#' @param B_matrices List with B1 and B0
#' @param m Maximum order
#' @param backend Character: "torch" (default) or "numpy"
#'
#' @return List with ATE, HOIF, and IIFF estimates for each order
#' @export
compute_hoif_estimators <- function(residuals, B_matrices, m = 5, backend = "torch") {
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

  # Compute U-statistics for each order j = 2 to m
  for (j in 2:m) {
    # Construct tensor list T_j^a
    # T_j^a = list(R^a, B^a, B^a, ..., B^a (j-1 times), r^a)
    T_j_1 <- list(R1)
    for (k in 1:(j-1)) {
      T_j_1[[k+1]] <- B1
    }
    T_j_1[[j+1]] <- r1

    T_j_0 <- list(R0)
    for (k in 1:(j-1)) {
      T_j_0[[k+1]] <- B0
    }
    T_j_0[[j+1]] <- r0

    # Construct expression E_j^a directly as nested list
    # For j tensors: [1,[1,2], [2,3], ..., [j-1,j],j]
    # This represents: R_i * B_{i,i1} * B_{i1,i2} * ... * r_{ij}
    E_j <- build_Ej(j)

    # Compute U-statistics using ustat function
    U1[j] <- (-1)^j * ustat(tensors = T_j_1, expression = E_j,
                            backend = backend, average = TRUE)
    U0[j] <- (-1)^j * ustat(tensors = T_j_0, expression = E_j,
                            backend = backend, average = TRUE)
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
#' @param pi Propensity scores
#' @param transform_method Character: "splines" or "fourier"
#' @param k Dimension of basis expansion
#' @param inverse_method Character: "direct", "nlshrink", or "corpcor"
#' @param m Maximum order for HOIF
#' @param sample_split Logical: whether to use sample splitting
#' @param K Number of folds for sample splitting (if used)
#' @param backend Character: "torch" (default) or "numpy"
#' @param seed Random seed for reproducibility (for sample splitting)
#' @param ... Additional arguments passed to transform_covariates
#'
#' @return List with ATE, HOIF, and IIFF estimates
#' @export
hoif_ate <- function(X, A, Y, mu1, mu0, pi,
                     transform_method = "splines",
                     k = 10,
                     inverse_method = "direct",
                     m = 5,
                     sample_split = FALSE,
                     K = 5,
                     backend = "torch",
                     seed = NULL,
                     ...) {

  n <- length(Y)

  # Step 1: Transform covariates (done on full data)
  Z <- transform_covariates(X, method = transform_method, k = k, ...)

  # Step 2: Compute residuals (done on full data)
  residuals <- compute_residuals(A, Y, mu1, mu0, pi)

  if (!sample_split) {
    # No sample splitting: standard procedure

    # Step 3: Compute inverse Gram matrices
    Omega <- compute_gram_inverse(Z, A, method = inverse_method)

    # Step 4: Compute basis matrices
    B_matrices <- compute_basis_matrix(Z, Omega$Omega1, Omega$Omega0)

    # Step 5: Compute HOIF estimators
    results <- compute_hoif_estimators(residuals, B_matrices, m, backend)

  } else {
    # Sample splitting (cross-fitting)

    # Set seed for fold creation if provided
    if (!is.null(seed)) {
      set.seed(seed)
    }

    # Create fold indices
    fold_indices <- sample(rep(1:K, length.out = n))

    # Storage for fold-specific estimates
    ATE_folds <- matrix(0, nrow = K, ncol = m - 1)
    HOIF1_folds <- matrix(0, nrow = K, ncol = m - 1)
    HOIF0_folds <- matrix(0, nrow = K, ncol = m - 1)
    IIFF1_folds <- matrix(0, nrow = K, ncol = m - 1)
    IIFF0_folds <- matrix(0, nrow = K, ncol = m - 1)

    for (j in 1:K) {
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

      results_j <- compute_hoif_estimators(residuals_j, B_matrices_j, m, backend)

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
