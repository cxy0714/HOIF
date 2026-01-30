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

# ==============================================================================
# Test 2: HOIF Estimation with Sample Splitting
# ==============================================================================
# This test evaluates the Higher-Order Influence Function (HOIF) estimator
# for the Average Treatment Effect (ATE) using sample splitting (cross-fitting).
#
# We simulate a partially linear data-generating process:
#   Y = A + X^T beta + noise
# where A is a binary treatment and X is a high-dimensional covariate vector.
#
# Nuisance outcome regressions E[Y | A=1, X] and E[Y | A=0, X] are estimated
# using Lasso, and the HOIF estimator is computed with and without sample splitting.
# We also verify reproducibility and backend consistency.
# ==============================================================================

cat("\n========================================\n")
cat("Test 2: HOIF main function \n")
cat("========================================\n")

# ------------------------------------------------------------------------------
# Step 1: Simulate synthetic data
# ------------------------------------------------------------------------------
set.seed(123)

n <- 1000 # sample size
p <- 50 # number of covariates

# Covariates
X <- matrix(rnorm(n * p), ncol = p)

# Binary treatment assignment (randomized)
A <- rbinom(n, 1, 0.5)

# True coefficient vector (normalized to unit length)
beta <- runif(p)
beta <- beta / sqrt(as.numeric(crossprod(beta)))

# Outcome: partially linear model
Y <- as.numeric(A + X %*% beta + rnorm(n, 0, 0.1))


# ------------------------------------------------------------------------------
# Step 2: Fit nuisance outcome regressions via Lasso
# ------------------------------------------------------------------------------
# We estimate:
#   mu1(X) = E[Y | A=1, X]
#   mu0(X) = E[Y | A=0, X]
# using separate Lasso regressions.

# Split data by treatment group
idx1 <- which(A == 1)
idx0 <- which(A == 0)

X1 <- X[idx1, , drop = FALSE]
Y1 <- Y[idx1]

X0 <- X[idx0, , drop = FALSE]
Y0 <- Y[idx0]

library(glmnet)

# Lasso regression for treated group
cv_fit1 <- cv.glmnet(X1, Y1, alpha = 1)
mu1_hat_all <- predict(cv_fit1, newx = X, s = "lambda.min")

# Lasso regression for control group
cv_fit0 <- cv.glmnet(X0, Y0, alpha = 1)
mu0_hat_all <- predict(cv_fit0, newx = X, s = "lambda.min")

mu1 <- as.vector(mu1_hat_all)
mu0 <- as.vector(mu0_hat_all)

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

m <- 7
n_folds <- 2

cat("Running HOIF with:\n")
cat("  n =", n, "\n")
cat("  p =", ncol(X), "\n")
cat("  order m =", m, "\n")
cat("  Sample split: TRUE (n_folds =", n_folds, ")\n")


# ------------------------------------------------------------------------------
# Step 4: Run HOIF with sample splitting (eHOIF)
# ------------------------------------------------------------------------------
results_split_1 <- hoif_ate(
  X, A, Y,
  mu1 = mu1,
  mu0 = mu0,
  pi = pi,
  transform_method = "none",
  m = m,
  sample_split = TRUE,
  n_folds = n_folds,
  seed = 123,
  backend = "torch"
)

cat("\nResults (with sample splitting):\n")
print(results_split_1)

# ==============================================================================
# Test Y: Double Check the eHOIF
# ==============================================================================

transform_method <- "none"
inverse_method <- "direct"
m <- 7
backend_1 <- "numpy"
backend_2 <- "torch"
n_folds <- 2

# Step 1: Transform covariates (done on full data)
Z <- transform_covariates(X, method = transform_method, basis_dim = p)

# Step 2: Compute residuals (done on full data)
residuals <- compute_residuals(A, Y, mu1, mu0, pi)


# Step X: Sample spliting
set.seed(123)
# Create fold indices
fold_indices <- sample(rep(1:n_folds, length.out = n))

I_1 <- which(fold_indices == 1)
I_2 <- which(fold_indices == 2)

Z_1 <- Z[I_1, , drop = FALSE]
A_1 <- A[I_1]
Z_2 <- Z[I_2, , drop = FALSE]
A_2 <- A[I_2]

residuals_1 <- list(
  R1 = residuals$R1[I_1],
  r1 = residuals$r1[I_1],
  R0 = residuals$R0[I_1],
  r0 = residuals$r0[I_1]
)
residuals_2 <- list(
  R1 = residuals$R1[I_2],
  r1 = residuals$r1[I_2],
  R0 = residuals$R0[I_2],
  r0 = residuals$r0[I_2]
)


# Step 3: Compute inverse Gram matrices
Omega_1 <- compute_gram_inverse(Z_2, A_2, method = inverse_method)
Omega_1_test <- (n / 2) * chol2inv(chol((t(A_2 * Z_2) %*% Z_2)))

norm(Omega_1_test - Omega_1$Omega1)

Omega_2_test <- (n / 2) * chol2inv(chol((t(A_1 * Z_1) %*% Z_1)))
Omega_2 <- compute_gram_inverse(Z_1, A_1, method = inverse_method)
norm(Omega_2_test - Omega_2$Omega1)


# Step 4: Compute basis matrices
B_matrices_1 <- compute_basis_matrix(Z_1, Omega_1$Omega1, Omega_2$Omega0)

B_matrices_1_test <- Z_1 %*% Omega_1_test %*% t(Z_1)
norm(B_matrices_1_test - B_matrices_1$B1)

B_matrices_2 <- compute_basis_matrix(Z_2, Omega_2$Omega1, Omega_1$Omega0)

B_matrices_2_test <- Z_2 %*% Omega_2_test %*% t(Z_2)
norm(B_matrices_2_test - B_matrices_2$B1)


# Step 5: Compute HOIF

results_1 <- compute_hoif_estimators(residuals_1, B_matrices_1, m, backend_2)
results_1$IIFF1

U2_1 <- ustat(list(residuals_1$R1, B_matrices_1_test, residuals_1$r1), "a,ab,b->", backend = "torch")
U3_1 <- ustat(list(residuals_1$R1, B_matrices_1_test, B_matrices_1_test, residuals_1$r1), "a,ab,bc,c->", backend = "torch")
U4_1 <- ustat(list(residuals_1$R1, B_matrices_1_test, B_matrices_1_test, B_matrices_1_test, residuals_1$r1), "a,ab,bc,cd,d->", backend = "torch")
U5_1 <- ustat(list(residuals_1$R1, B_matrices_1_test, B_matrices_1_test, B_matrices_1_test, B_matrices_1_test, residuals_1$r1), "a,ab,bc,cd,de,e->", backend = "torch")
U6_1 <- ustat(list(residuals_1$R1, B_matrices_1_test, B_matrices_1_test, B_matrices_1_test, B_matrices_1_test, B_matrices_1_test, residuals_1$r1), "a,ab,bc,cd,de,ef,f->", backend = "torch")

IF2_1 <- U2_1
IF3_1 <- -U3_1 + U2_1
IF4_1 <- U4_1 - 2 * U3_1 + U2_1
IF5_1 <- -U5_1 + 3 * U4_1 - 3 * U3_1 + U2_1
c(IF2_1, IF3_1, IF4_1, IF5_1)
results_1$IIFF1

results_2 <- compute_hoif_estimators(residuals_2, B_matrices_2, m, backend_2)
results_2$IIFF1

U2_2 <- ustat(list(residuals_2$R1, B_matrices_2_test, residuals_2$r1), "a,ab,b->", backend = "torch")
U3_2 <- ustat(list(residuals_2$R1, B_matrices_2_test, B_matrices_2_test, residuals_2$r1), "a,ab,bc,c->", backend = "torch")
U4_2 <- ustat(list(residuals_2$R1, B_matrices_2_test, B_matrices_2_test, B_matrices_2_test, residuals_2$r1), "a,ab,bc,cd,d->", backend = "torch")
U5_2 <- ustat(list(residuals_2$R1, B_matrices_2_test, B_matrices_2_test, B_matrices_2_test, B_matrices_2_test, residuals_2$r1), "a,ab,bc,cd,de,e->", backend = "torch")
U6_2 <- ustat(list(residuals_2$R1, B_matrices_2_test, B_matrices_2_test, B_matrices_2_test, B_matrices_2_test, B_matrices_2_test, residuals_2$r1), "a,ab,bc,cd,de,ef,f->", backend = "torch")
U7_2 <- ustat(list(residuals_2$R1, B_matrices_2_test, B_matrices_2_test, B_matrices_2_test, B_matrices_2_test, B_matrices_2_test, residuals_2$r1), "a,ab,bc,cd,de,ef,f->", backend = "torch")

IF2_2 <- +U2_2
IF3_2 <- -U3_2 + U2_2
IF4_2 <- +U4_2 - 2 * U3_2 + U2_2
IF5_2 <- -U5_2 + 3 * U4_2 - 3 * U3_2 + U2_2
IF6_2 <- +U6_2 - 4 * U5_2 + 6 * U4_2 - 4 * U3_2 + U2_2
IF7_2 <- -U7_2 - 6 * U6_2 - 15 * U5_2 + 6 * U4_2 - 4 * U3_2 + U2_2

c(IF2_2, IF3_2, IF4_2, IF5_2, IF6_2, IF7_2)
results_2$IIFF1


# Step 5: Compute HOIF
(results_2$IIFF1 + results_1$IIFF1) / 2
results_split_1$IIFF1

# ------------------------------------------------------------------------------
# Step 5: Reproducibility check
# ------------------------------------------------------------------------------
# Running the estimator again with the same seed should produce identical results.
results_split2 <- hoif_ate(
  X, A, Y,
  mu1 = mu1,
  mu0 = mu0,
  pi = pi,
  transform_method = "none",
  m = m,
  sample_split = TRUE,
  n_folds = n_folds,
  seed = 123,
  backend = "torch"
)

cat("\nReproducibility check:\n")
cat("  Same ATE:", all.equal(results_split_1$ATE, results_split2$ATE), "\n")

stopifnot(all.equal(results_split_1$ATE, results_split2$ATE))

cat("\n✓ Sample splitting is reproducible!\n")


# ------------------------------------------------------------------------------
# Step 6: Run HOIF without sample splitting (sHOIF)
# ------------------------------------------------------------------------------

results_non_sample_splitting <- hoif_ate(
  X, A, Y,
  mu1 = mu1,
  mu0 = mu0,
  pi = pi,
  transform_method = "none",
  m = m,
  sample_split = 0,
  backend = "torch"
)

cat("\nResults without sample splitting:\n")
print(results_non_sample_splitting)
plot(results_non_sample_splitting)

cat("\n✓ HOIF without sample splitting works!\n")


# ------------------------------------------------------------------------------
# Step 7: Backend consistency check (Torch vs Pure R)
# ------------------------------------------------------------------------------
results_non_sample_splitting_pure_R <- hoif_ate(
  X, A, Y,
  mu1 = mu1,
  mu0 = mu0,
  pi = pi,
  transform_method = "none",
  m = 6,
  sample_split = 0,
  backend = "torch",
  pure_R_code = TRUE
)

cat("\nResults with pure R backend (no sample splitting):\n")
print(results_non_sample_splitting_pure_R)
plot(results_non_sample_splitting_pure_R)



# ==============================================================================
# Test 3: Call all intermediate functions include ustat()
# ==============================================================================
cat("\n========================================\n")
cat("Test 3: Call all intermediate functions include ustat()\n")
cat("========================================\n")


transform_method <- "none"
inverse_method <- "direct"
m <- 6
backend_1 <- "numpy"
backend_2 <- "torch"

# Step 1: Transform covariates (done on full data)
Z <- transform_covariates(X, method = transform_method, basis_dim = p)

# Step 2: Compute residuals (done on full data)
residuals <- compute_residuals(A, Y, mu1, mu0, pi)

# Step 3: Compute inverse Gram matrices
Omega <- compute_gram_inverse(Z, A, method = inverse_method)

# Step 4: Compute basis matrices
B_matrices <- compute_basis_matrix(Z, A, Omega$Omega1, Omega$Omega0)

# Step 5: Compute HOIF estimators
results_1 <- compute_hoif_estimators(residuals, B_matrices, m, backend_1)
results_2 <- compute_hoif_estimators(residuals, B_matrices, m, backend_2)

# When using pure_R_code, whether m is, always compute 6-th order
results_3 <- compute_hoif_estimators(residuals, B_matrices, m, backend = "torch", pure_R_code = TRUE)

stopifnot(all.equal(results_1, results_2, tolerance = 1e-6))
stopifnot(all.equal(results_2, results_3, tolerance = 1e-6))
stopifnot(all.equal(results_1, results_3, tolerance = 1e-6))

# Test ustat
U2 <- ustat(list(residuals$R1, B_matrices$B1, residuals$r1), "a,ab,b->", backend = "torch")
U3 <- ustat(list(residuals$R1, B_matrices$B1, B_matrices$B1, residuals$r1), "a,ab,bc,c->", backend = "torch")
U4 <- ustat(list(residuals$R1, B_matrices$B1, B_matrices$B1, B_matrices$B1, residuals$r1), "a,ab,bc,cd,d->", backend = "torch")
U5 <- ustat(list(residuals$R1, B_matrices$B1, B_matrices$B1, B_matrices$B1, B_matrices$B1, residuals$r1), "a,ab,bc,cd,de,e->", backend = "torch")
U6 <- ustat(list(residuals$R1, B_matrices$B1, B_matrices$B1, B_matrices$B1, B_matrices$B1, B_matrices$B1, residuals$r1), "a,ab,bc,cd,de,ef,f->", backend = "torch")

cat("\n✓ ustat works with both formats!\n")

R_u <- calculate_u_statistics_six(Vector_1 = residuals$R1, Vector_2 = residuals$r1, A1 = B_matrices$B1, A2 = B_matrices$B1, A3 = B_matrices$B1, A4 = B_matrices$B1, A5 = B_matrices$B1)
str(R_u)
str(c(U2, U3, U4, U5, U6))


# ==============================================================================
# Test X: Expression Conversion
# ==============================================================================
cat("\n========================================\n")
cat("Test X: Expression Conversion\n")
cat("========================================\n")

# Test the helper function
expr1 <- HOIF:::build_Ej(3)
result1 <- HOIF:::expr_list_to_einstein(expr1)
print(str(expr1))
cat("Output:", result1, "\n")
cat("Expected: a,ab,bc,c->\n")
stopifnot(result1 == "a,ab,bc,c->")
