source("R/hoif_ustat.R")

# package
library("SMUT") # for efficient matrix multiplication via eigenMapMatMult()

HOIF <- function(A, Y, X, a_est, b_est_1, b_est_0, order, k,
                 Z_method = "spline",
                 Omega_method = "robust",
                 is_split = FALSE,
                 is_bootstrap = FALSE,
                 bootstrap_seed = 42,
                 bootstrap_number = 1000) {

  # 参数验证
  validate_inputs(A, Y, X, a_est, b_est_1,b_est_0, order, k, Z_method, Omega_method,
                  is_split, is_bootstrap, bootstrap_seed, bootstrap_number)

  n <- length(A)

  Z_k <- Z_transform(k, X, Z_method)
  e_A <- Epsilon_A(A, a_est = a_est)
  e_Y <- Epsilon_Y(Y, A, b_est_1 = b_est_1, b_est_0 = b_est_0)

  if (is_split) {
    split_indices <- sample(1:n, n/2, replace = FALSE)
    basis <- Basis_Omega_estimation(Z_k = Z_k, A = A,
                              Omega_method = Omega_method,
                              is_split = is_split,
                              split_indices = split_indices)
  } else {
    basis <- Basis_Omega_estimation(Z_k = Z_k, A = A,
                              Omega_method = Omega_method)
  }
  basis <-

  IIFF <- IIFF_estimation(Omega, Z_k, e_A, e_Y, order, is_split)

  if (is_bootstrap) {
    set.seed(bootstrap_seed)

    weights <- matrix(
      sample(1:n, n * bootstrap_number, replace = TRUE),
      nrow = n,
      ncol = bootstrap_number
    )

    IIFF_bootstrap <- list()
    for (i in 1:bootstrap_number) {
      IIFF_bootstrap[[i]] <- IIFF_estimation(Omega, Z_k, e_A, e_Y,
                                             weights = weights[, i], order)
    }
    var_bootstrap <- IIFF_var_estimation(IIFF_bootstrap)
  }


  result <- format_results(IIFF, a, is_bootstrap,
                           if(is_bootstrap) var_bootstrap else NULL)

  return(result)
}


format_results <- function(IIFF, a, is_bootstrap, var_est = NULL) {

  if (length(a) == 2 && all(a == c(0, 1))) {
    result <- list(
      ATE_bias = IIFF$a_1 - IIFF$a_0,
      Y_1_bias = IIFF$a_1,
      Y_0_bias = IIFF$a_0
    )

    if (is_bootstrap) {
      result$var_est <- var_est
      result$se_est <- sqrt(var_est)
    }
  } else {

    result <- list(
      Y_bias = IIFF,
      a_value = a
    )
  }

  class(result) <- "HOIF"
  return(result)
}


Z_transform <- function(k, X, Z_method) {
  return( Z_k = X)
}

Basis_Omega_estimation <- function(Z_k, A, a, Omega_method,
                             is_split = FALSE, split_indices = NULL) {

  if (is_split){
    Z_k[split_indices]
  }
}

Epsilon_A <- function(A, a_est) {
    results <- list()
    results$a_1 <-  1 - A/a_est
    results$a_0 <-  1 - (1-A)/(1-a_est)
  return(results)
}

Epsilon_Y <- function(Y, A, b_est_1, b_est_0) {
    results <- list()
    results$a_1 <-  A(Y - b_est_1)
    results$a_0 <-  (1-A)(Y - b_est_0)
    return(results)
}

IIFF_estimation <- function(Omega, Z_k, e_A, e_Y, order, weights = NULL) {

}

IIFF_estimation_single <- function(Basis, e_A, e_Y, order, weights = NULL) {

}
IIFF_var_estimation <- function(IIFF_bootstrap) {
  # 计算bootstrap方差
  # 实现
}

################################################################################
# Main function using "all" U-statistics definition! i.e. (i_1 \ne i_2 \ne ... \ne i_m)
# Computes and return 4's estimators from 2-th HOIF to 6-th HOIF estimator. Higher order code is coming soon!
# Vector_1: a n-dimensional vector containing the treatment residuals (like (Aa -1) )
# Vector_2: a n-dimensional vector containing the outcome residuals (like (y - b),here no A-weighted)
# weight : a n-dimensional vector containing the treatment for Weight-ed Gram matrix ( like A or (1-A) )
# basis: a n*p matrix containing the basis transformations of the confounder.
# order: number of estimators, default value is 6. only can be 2 or 3 or 4 or 5 or 6.
# Split: a logistic variable, default value is 0.
#        Split== 1 means using sample Splitto compute the Omega_hat (eHOIF),
#        Split== 0 means using whole sample to compute the Omega_hat (sHOIF).
################################################################################
compute_HOIF_general_all_U <- function(Vector_1, Vector_2, weight, basis, order = 6, Split= 0) {
  n <- length(Vector_1)

  if ( Split== 0) {
    basis_weight <- basis * weight
    t_basis_weight <- t(basis_weight)
    Mat <- eigenMapMatMult(t_basis_weight, basis)
    L <- chol(Mat)
    Omega_mat <- chol2inv(L) * n


    Ker <- eigenMapMatMult(eigenMapMatMult(basis, Omega_mat), t_basis_weight)

    U_list <- list()
    results <- calculate_u_statistics_six(Vector_1 = Vector_1, Vector_2 = Vector_2, A1 = Ker, A2 = Ker, A3 = Ker, A4 = Ker, A5 = Ker)

    for (i in 2:(order)) {
      U_list[[paste0("U_", i)]] <- (-1)^i * results[[i-1]]
    }


    BD_list <- Buildingblock_sum_HOIF(U_list, order = order)

    return_list <- c(BD_list, U_list)
    return_list <- combine_results(sHOIF = return_list)
  } else if ( Split== 1) {


    # Splitindices into two halves: a and b
    indices <- seq_len(n)
    set.seed(123)  # For reproducibility
    indices_a <- sample(indices, round(n / 2))
    indices_b <- setdiff(indices, indices_a)

    # Splitbasis, weight, Vector_1, Vector_2
    basis_a <- basis[indices_a, ]
    basis_b <- basis[indices_b, ]
    weight_a <- weight[indices_a]
    weight_b <- weight[indices_b]
    Vector_1_b <- Vector_1[indices_b]
    Vector_2_b <- Vector_2[indices_b]

    # Compute for group a
    basis_weight_a <- basis_a * weight_a
    t_basis_weight_a <- t(basis_weight_a)
    Mat_a <- eigenMapMatMult(t_basis_weight_a, basis_a)
    L_a <- chol(Mat_a)
    Omega_mat_a <- chol2inv(L_a) * (n / 2)

    # Compute Ker using Omega_mat_a and basis_b
    Ker_ab <- eigenMapMatMult(eigenMapMatMult(basis_b, Omega_mat_a), t(basis_weight_a))

    # Compute results using Ker_ab and Vector_1_b, Vector_2_b
    results_ab <- calculate_u_statistics_six(Vector_1 = Vector_1_b, Vector_2 = Vector_2_b,
                                             A1 = Ker_ab, A2 = Ker_ab, A3 = Ker_ab, A4 = Ker_ab,A5 = Ker_ab)

    U_list_ab <- list()
    for (i in 2:order) {
      U_list_ab[[paste0("U_", i)]] <- (-1)^i * results_ab[[i - 1]]
    }

    BD_list_ab <- Buildingblock_sum_HOIF(U_list_ab, order = order)
    return_list_ab <- c(BD_list_ab, U_list_ab)

    # Combine results for SH-OIF (Splita -> b)
    return_list_ab <- combine_results(eHOIF = return_list_ab)

    # Repeat the process, swapping a and b
    basis_weight_b <- basis_b * weight_b
    t_basis_weight_b <- t(basis_weight_b)
    Mat_b <- eigenMapMatMult(t_basis_weight_b, basis_b)
    L_b <- chol(Mat_b)
    Omega_mat_b <- chol2inv(L_b) * (n / 2)

    Ker_ba <- eigenMapMatMult(eigenMapMatMult(basis_a, Omega_mat_b), t(basis_weight_b))

    results_ba <- calculate_u_statistics_six(Vector_1 = Vector_1[indices_a], Vector_2 = Vector_2[indices_a],
                                             A1 = Ker_ba, A2 = Ker_ba, A3 = Ker_ba, A4 = Ker_ba, A5 = Ker_ba)

    U_list_ba <- list()
    for (i in 2:order) {
      U_list_ba[[paste0("U_", i)]] <- (-1)^i * results_ba[[i - 1]]
    }

    BD_list_ba <- Buildingblock_sum_HOIF(U_list_ba, order = order)
    return_list_ba <- c(BD_list_ba, U_list_ba)

    return_list_ba <- combine_results(eHOIF = return_list_ba)

    # Average the two return lists
    return_list <- add_lists(return_list_ab, return_list_ba)
    return_list <- lapply(return_list, function(x) x / 2)

  }

  return(return_list)
}


################################################################################
# Main function using "part" U-statistics definition i.e. (i_1 < i_2 < ... < i_m)
# Computes and return any order>2 HOIF estimator,(order-1)'s estimators from 2-th HOIF to order-th HOIF estimator .
# Vector_1: a n-dimensional vector containing the treatment residuals (like (Aa -1) )
# Vector_2: a n-dimensional vector containing the outcome residuals (like (y - b) ,but not (A(y-b))! Here no A-weighted)
# weight : a n-dimensional vector containing the treatment for Weight-ed Gram matrix ( like A or (1-A) )
# basis: a n*p matrix containing the basis transformations of the confounder.
# order: number of estimators.
# Split: a logistic variable, default value is 0.
#        Split== 1 means using sample Splitto compute the Omega_hat (eHOIF),
#        Split== 0 means using whole sample to compute the Omega_hat (sHOIF).
################################################################################
compute_HOIF_general_part_U <- function(Vector_1, Vector_2, weight, basis, order, Split= 0) {
  n <- length(Vector_1)

  if (Split== 0) {
    # Original implementation for Split== 0
    basis_weight <- basis * weight
    t_basis_weight <- t(basis_weight)
    Mat <- eigenMapMatMult(t_basis_weight, basis)
    L <- chol(Mat)
    Omega_mat <- chol2inv(L) * n

    Ker <- eigenMapMatMult(eigenMapMatMult(basis, Omega_mat), t_basis_weight)
    Ker_up <- Ker * upper.tri(Ker)
    t_vector_1 <- t(Vector_1)
    C <- diag(n)

    U_list <- list()

    for (i in 2:order) {
      C <- eigenMapMatMult(C, Ker_up)
      U_list[[paste0("U_", i)]] <- (-1)^i * as.numeric(eigenMapMatMult(eigenMapMatMult(t_vector_1, C), Vector_2)) / choose(n, i)
    }

    BD_list <- Buildingblock_sum_HOIF(U_list, order)
    return_list <- c(BD_list, U_list)
    return_list <- combine_results(sHOIF = return_list)

  } else if (Split== 1) {
    # Sample-Splitimplementation for Split== 1
    # Splitindices into two halves: a and b
    idx_a <- sample(1:n, size = round(n / 2), replace = FALSE)
    idx_b <- setdiff(1:n, idx_a)

    # Splitbasis, weight, Vector_1, and Vector_2
    basis_a <- basis[idx_a, , drop = FALSE]
    basis_b <- basis[idx_b, , drop = FALSE]
    weight_a <- weight[idx_a]
    weight_b <- weight[idx_b]
    Vector_1_a <- Vector_1[idx_a]
    Vector_2_a <- Vector_2[idx_a]
    Vector_1_b <- Vector_1[idx_b]
    Vector_2_b <- Vector_2[idx_b]

    # Compute Mat and Omega_mat using subset a
    basis_weight_a <- basis_a * weight_a
    t_basis_weight_a <- t(basis_weight_a)
    Mat_a <- eigenMapMatMult(t_basis_weight_a, basis_a)
    L_a <- chol(Mat_a)
    Omega_mat_a <- chol2inv(L_a) * (n / 2)

    # Compute Ker using basis_b and Omega_mat_a
    t_basis_weight_b <- t(basis_b * weight_b)
    Ker_ab <- eigenMapMatMult(eigenMapMatMult(basis_b, Omega_mat_a), t_basis_weight_b)
    Ker_up_ab <- Ker_ab * upper.tri(Ker_ab)

    # Calculate U_list for subset b
    t_vector_1_b <- t(Vector_1_b)
    C_ab <- diag(length(idx_b))
    U_list_ab <- list()

    for (i in 2:order) {
      C_ab <- eigenMapMatMult(C_ab, Ker_up_ab)
      U_list_ab[[paste0("U_", i)]] <- (-1)^i * as.numeric(eigenMapMatMult(eigenMapMatMult(t_vector_1_b, C_ab), Vector_2_b)) / choose((n / 2), i)
    }

    # Compute Mat and Omega_mat using subset b
    basis_weight_b <- basis_b * weight_b
    t_basis_weight_b <- t(basis_weight_b)
    Mat_b <- eigenMapMatMult(t_basis_weight_b, basis_b)
    L_b <- chol(Mat_b)
    Omega_mat_b <- chol2inv(L_b) * (n / 2)

    # Compute Ker using basis_a and Omega_mat_b
    t_basis_weight_a <- t(basis_a * weight_a)
    Ker_ba <- eigenMapMatMult(eigenMapMatMult(basis_a, Omega_mat_b), t_basis_weight_a)
    Ker_up_ba <- Ker_ba * upper.tri(Ker_ba)

    # Calculate U_list for subset a
    t_vector_1_a <- t(Vector_1_a)
    C_ba <- diag(length(idx_a))
    U_list_ba <- list()

    for (i in 2:order) {
      C_ba <- eigenMapMatMult(C_ba, Ker_up_ba)
      U_list_ba[[paste0("U_", i)]] <- (-1)^i * as.numeric(eigenMapMatMult(eigenMapMatMult(t_vector_1_a, C_ba), Vector_2_a)) / choose((n / 2), i)
    }

    # Combine results from both splits and average
    BD_list_ab <- Buildingblock_sum_HOIF(U_list_ab, order)
    BD_list_ba <- Buildingblock_sum_HOIF(U_list_ba, order)

    return_list_ab <- c(BD_list_ab, U_list_ab)
    return_list_ab <- combine_results(eHOIF = return_list_ab)

    return_list_ba <- c(BD_list_ba, U_list_ba)
    return_list_ba <- combine_results(eHOIF = return_list_ba)

    # Average the two return lists
    return_list <- add_lists(return_list_ab, return_list_ba)
    return_list <- lapply(return_list, function(x) x / 2)
  }

  return(return_list)
}


################################################################################
# Helper funciton for compute_HOIF_general_part_U()
# Computes and returns (order)-th HOIF estimator and IF from buildingblock U_list
################################################################################
Buildingblock_sum_HOIF <- function(U_list, order, fplugin = 0) {
  BD_list <- list()

  for (j in 2:(order)) {
    BD_value <- 0

    for (i in 0:(j - 2)) {

      combination <- choose(j - 2, i)

      U_index <- j - i

      BD_value <- BD_value + combination * U_list[[paste0("U_", U_index)]]
    }


    BD_list[[paste0("IF_", j)]] <- BD_value
  }


  sum_list <- list()


  for (i in 2:order) {
    sum_list[[paste0("HOIF", i)]] <- sum(unlist(BD_list[1:(i - 1)])) + fplugin
  }
  return_list <- c(sum_list, BD_list)
  return(return_list)
}
################################################################################
# Main body function in compute_HOIF_general_all_U()
# Computes and return 4 building_block estimators from 2-th HOIF to 6-th U-statistics.
# \bbU_{n, m}^{\mathsf{all}} \biggl[ f(O_{i_1}, \ldots, O_{i_m}) \biggr] = \frac{1}{\binom{n}{m} \factorial{m} } \sum_{i_1 \ne i_2 \ne \ldots \ne i_m} f(O_{i_1}, \ldots, O_{i_m})
################################################################################

################################################################################
# Helper funciton for compute_HOIF_general_part_U() and compute_HOIF_general_all_U()
# subtract, add, combine function for list
################################################################################
# Helper function to create no-diagonal matrix
no_diag <- function(mat) {
  diag(mat) <- 0
  return(mat)
}

# Helper function for Hadamard (point-wise) product
hadamard <- function(mat1, mat2) {
  return(mat1 * mat2)
}

# Helper function for Diag-Column-Sum
diag_col_sum <- function(mat) {
  diag(colSums(mat))
}
# Helper function for Diag-Column-Sum
diag_Mat <- function(mat) {
  diag(diag(mat))
}

subtract_lists <- function(list1, list2) {
  if (length(list1) != length(list2)) {
    stop("Lists must be of the same length")
  }
  mapply(function(x, y) {
    x - y
  }, list1, list2, SIMPLIFY = FALSE)
}
add_lists <- function(list1, list2) {
  if (length(list1) != length(list2)) {
    stop("Lists must be of the same length")
  }
  mapply(function(x, y) {
    x + y
  }, list1, list2, SIMPLIFY = FALSE)
}
combine_results <- function(...) {

  input_lists <- list(...)
  list_names <- names(input_lists)


  if (is.null(list_names) || any(list_names == "")) {
    stop("All input lists must be named")
  }


  processed_lists <- lapply(names(input_lists), function(list_name) {
    current_list <- input_lists[[list_name]]

    if (!is.list(current_list)) {
      current_list <- list(current_list)
    }


    result <- lapply(names(current_list), function(name) {
      current_list[[name]]
    })
    names(result) <- paste0(list_name, "_", names(current_list))

    return(result)
  })


  result <- do.call(c, processed_lists)
  return(result)
}

################## Example ######################
# In this example, we estimate the bias of DML/AIPW for the potential outcomes Y(a = 1) and Y(a = 0).
# The outcome model assumes a linear relationship: Y = A + beta * X.
# The propensity score model follows a logistic regression: A = psi(alpha * X).
#
# - compute_HOIF_general_all_U(): Provides the exact formula for HOIF, but only up to the 5th order.
# - compute_HOIF_general_part_U(): Computes an approximate formula for HOIF, extendable to any order.
#       - in their return results, "_HOIF_" are the used estimator, "_IF_" or "_U_" are just middle terms.
# Note: The basis function of X is the identity in this example. For practical applications,
#       it is recommended to apply transformations (e.g., B-splines, Fourier basis, etc.) to X.
################################################

set.seed(123)


n <- 1000
p <- 10

mu_x <- 0
alpha <- rnorm(p, mean = 0, sd = 0.1)
beta <- rnorm(p, mean = 0, sd = 0.1)

X <- matrix(rnorm(n * p, mu_x, sd = 1), nrow = n, ncol = p)

logit_prob <- X %*% alpha
prob_A <- 1 / (1 + exp(-logit_prob))
A <- as.numeric(rbinom(n, size = 1, prob = prob_A))


Y <- as.numeric(A + X %*% beta + rnorm(n, mean = 0, sd = 0.2))


propensity_model <- glm(A ~ X, family = binomial(link = "logit"))
propensity_scores <- predict(propensity_model, type = "response")
propensity_scores <- as.numeric(propensity_scores)


Y_model <- lm(Y ~ A + X)


Y_pred_1 <- predict(Y_model, newdata = data.frame(A = 1, X = X))
Y_pred_0 <- predict(Y_model, newdata = data.frame(A = 0, X = X))


summary(propensity_model)
summary(Y_model)

epsilon_A_1 <- A / propensity_scores - 1
epsilon_A_0 <- (1 - A) / (1 - propensity_scores) - 1
epsilon_Y_1 <- Y - Y_pred_1
epsilon_Y_0 <- Y - Y_pred_0

m_part <- 10
m_all <- 6

hoif_1_all_shoif <- compute_HOIF_general_all_U(
  Vector_1 = epsilon_A_1,
  Vector_2 = epsilon_Y_1,
  weight = A,
  basis = X,
  order = m_all,
  Split= 0
)
hoif_1_part_shoif <- compute_HOIF_general_part_U(
  Vector_1 = epsilon_A_1,
  Vector_2 = epsilon_Y_1,
  weight = A,
  basis = X,
  order = m_part,
  Split= 0
)

hoif_1_all_ehoif <- compute_HOIF_general_all_U(
  Vector_1 = epsilon_A_1,
  Vector_2 = epsilon_Y_1,
  weight = A,
  basis = X,
  order = m_all,
  Split= 1
)
hoif_1_part_ehoif <- compute_HOIF_general_part_U(
  Vector_1 = epsilon_A_1,
  Vector_2 = epsilon_Y_1,
  weight = A,
  basis = X,
  order = m_part,
  Split= 1
)


hoif_0_all_shoif <- compute_HOIF_general_all_U(
  Vector_1 = epsilon_A_0,
  Vector_2 = epsilon_Y_0,
  weight = 1 - A,
  basis = X,
  order = m_all,
  Split= 0
)
hoif_0_part_shoif <- compute_HOIF_general_part_U(
  Vector_1 = epsilon_A_0,
  Vector_2 = epsilon_Y_0,
  weight = 1 - A,
  basis = X,
  order = m_part,
  Split= 0
)

hoif_0_all_ehoif <- compute_HOIF_general_all_U(
  Vector_1 = epsilon_A_0,
  Vector_2 = epsilon_Y_0,
  weight = 1 - A,
  basis = X,
  order = m_all,
  Split= 1
)

hoif_0_part_ehoif <- compute_HOIF_general_part_U(
  Vector_1 = epsilon_A_0,
  Vector_2 = epsilon_Y_0,
  weight = 1 - A,
  basis = X,
  order = m_part,
  Split= 1
)
