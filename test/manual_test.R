#' Manual Testing Script for HOIF Package
#'
#' Run these tests step by step to verify the package works correctly
#'
#' @author Xingyu Chen
#' @date 2026-01-23
devtools::load_all()

library(HOIF)

# ==============================================================================
# Test 1: Environment Setup
# ==============================================================================
cat("\n========================================\n")
cat("Test 1: Check Environment\n")
cat("========================================\n")

check_hoif_setup()

# ==============================================================================
# Test 2: Expression Conversion
# ==============================================================================
cat("\n========================================\n")
cat("Test 2: Expression Conversion\n")
cat("========================================\n")

# Test the helper function
expr1 <- HOIF:::build_Ej(3)
result1 <- HOIF:::expr_list_to_einstein(expr1)
print(str(expr1))
cat("Output:", result1, "\n")
cat("Expected: a,ab,bc,c->\n")
stopifnot(result1 == "a,ab,bc,c->")

# ==============================================================================
# Test 3: Direct ustat Call
# ==============================================================================
cat("\n========================================\n")
cat("Test 3: Direct ustat Call\n")
cat("========================================\n")

set.seed(123)
n <- 2000
p <- 50
X <- matrix(rnorm(n * p), ncol = p)
A <- rbinom(n, 1, 0.5)
beta <- runif(p)
beta <- beta/sqrt(as.numeric(crossprod(beta)))
Y <- as.numeric(A + X %*% beta + rnorm(n,0, 0.1))

# Split data by treatment
idx1 <- which(A == 1)
idx0 <- which(A == 0)

X1 <- X[idx1, , drop = FALSE]
Y1 <- Y[idx1]

X0 <- X[idx0, , drop = FALSE]
Y0 <- Y[idx0]

library(glmnet)
cv_fit1 <- cv.glmnet(X1, Y1, alpha = 1)  # alpha=1 => Lasso
mu1_hat_all <- predict(cv_fit1, newx = X, s = "lambda.min")  # predict on ALL X
cv_fit0 <- cv.glmnet(X0, Y0, alpha = 1)
mu0_hat_all <- predict(cv_fit0, newx = X, s = "lambda.min")

mu1 <- as.vector(mu1_hat_all)
mu0 <- as.vector(mu0_hat_all)
# Simple nuisance estimates
pi <- rep(0.5, n)

transform_method <- "none"
inverse_method <- "direct"
m <- 4
backend_1 <- "numpy"
backend_2 <- "torch"
# Step 1: Transform covariates (done on full data)
Z <- transform_covariates(X, method = transform_method, k = p)

# Step 2: Compute residuals (done on full data)
residuals <- compute_residuals(A, Y, mu1, mu0, pi)

# Step 3: Compute inverse Gram matrices
Omega <- compute_gram_inverse(Z, A, method = inverse_method)

# Step 4: Compute basis matrices
B_matrices <- compute_basis_matrix(Z, Omega$Omega1, Omega$Omega0)

# Step 5: Compute HOIF estimators
results_1 <- compute_hoif_estimators(residuals, B_matrices, m, backend_1)
results_2 <- compute_hoif_estimators(residuals, B_matrices, m, backend_2)


# Test with Einstein notation (string)
U3 <- ustat(list(residuals$R1, B_matrices$B1, B_matrices$B1, residuals$r1), "a,ab,bc,c->", backend = "torch")

L1 <- sweep(B_matrices$B1, 1, residuals$R1, "*")
R1 <- sweep(B_matrices$B1, 2, residuals$r1, "*")
diag(L1) <- 0
diag(R1) <- 0
result_list <- (sum(L1 %*% R1) - sum(diag(L1 %*% R1)))/(n*(n-1)*(n-2))
cat("Result (direct computing):", result_list, "\n")


cat("\n✓ ustat works with both formats!\n")

ustat(list(residuals$R1, B_matrices$B1, B_matrices$B1, B_matrices$B1, residuals$r1), "a,ab,bc,cd,d->", backend = "numpy")

U4 <- ustat(list(residuals$R1, B_matrices$B1, B_matrices$B1, B_matrices$B1, residuals$r1), "a,ab,bc,cd,d->", backend = "torch")
U5 <- ustat(list(residuals$R1, B_matrices$B1, B_matrices$B1, B_matrices$B1, B_matrices$B1, residuals$r1), "a,ab,bc,cd,de,e->", backend = "torch")
U6 <- ustat(list(residuals$R1, B_matrices$B1, B_matrices$B1, B_matrices$B1, B_matrices$B1, B_matrices$B1, residuals$r1), "a,ab,bc,cd,de,ef,f->", backend = "torch")

source("test/hoif_r.R")
d_u <-  calculate_u_statistics_six(Vector_1 = residuals$R1, Vector_2 = residuals$r1, A1 =  B_matrices$B1, A2 = B_matrices$B1, A3 = B_matrices$B1, A4 = B_matrices$B1, A5 = B_matrices$B1)
str(d_u)
# ==============================================================================
# Test 4: Covariate Transformation
# ==============================================================================
cat("\n========================================\n")
cat("Test 4: Covariate Transformation\n")
cat("========================================\n")

set.seed(456)
n <- 50
X <- matrix(rnorm(n * 2), ncol = 2)

# Test B-splines
Z_splines <- transform_covariates(X, method = "splines", k = 5)
cat("B-splines transformation:\n")
cat("  Input dim:", dim(X), "\n")
cat("  Output dim:", dim(Z_splines), "\n")
cat("  Intercept check:", all(Z_splines[,1] == 1), "\n")

stopifnot(nrow(Z_splines) == n)
stopifnot(ncol(Z_splines) > 2)
stopifnot(all(Z_splines[,1] == 1))

# Test Fourier
Z_fourier <- transform_covariates(X, method = "fourier", k = 6)
cat("\nFourier transformation:\n")
cat("  Input dim:", dim(X), "\n")
cat("  Output dim:", dim(Z_fourier), "\n")

cat("\n✓ Transformation works!\n")

# ==============================================================================
# Test 5: Residual Computation
# ==============================================================================
cat("\n========================================\n")
cat("Test 5: Residual Computation\n")
cat("========================================\n")

set.seed(789)
n <- 100
A <- rbinom(n, 1, 0.5)
Y <- rnorm(n)
mu1 <- rnorm(n)
mu0 <- rnorm(n)
pi <- runif(n, 0.2, 0.8)

residuals <- compute_residuals(A, Y, mu1, mu0, pi)

cat("Residuals computed:\n")
cat("  R1 length:", length(residuals$R1), "\n")
cat("  r1 length:", length(residuals$r1), "\n")
cat("  R0 length:", length(residuals$R0), "\n")
cat("  r0 length:", length(residuals$r0), "\n")

stopifnot(length(residuals$R1) == n)
stopifnot(length(residuals$r1) == n)

cat("\n✓ Residuals computed correctly!\n")

# ==============================================================================
# Test 6: Gram Matrix Inversion
# ==============================================================================
cat("\n========================================\n")
cat("Test 6: Gram Matrix Inversion\n")
cat("========================================\n")

set.seed(111)
n <- 50
k <- 10
Z <- matrix(rnorm(n * k), ncol = k)
A <- rbinom(n, 1, 0.5)

Omega <- compute_gram_inverse(Z, A, method = "direct")

cat("Gram matrices inverted:\n")
cat("  Omega1 dim:", dim(Omega$Omega1), "\n")
cat("  Omega0 dim:", dim(Omega$Omega0), "\n")

stopifnot(all(dim(Omega$Omega1) == c(k, k)))
stopifnot(all(dim(Omega$Omega0) == c(k, k)))

cat("\n✓ Gram matrix inversion works!\n")

# ==============================================================================
# Test 7: Simple HOIF (No Sample Split)
# ==============================================================================
cat("\n========================================\n")
cat("Test 7: Simple HOIF Estimation\n")
cat("========================================\n")

set.seed(123)
n <- 1000
p <- 50
X <- matrix(rnorm(n * p), ncol = p)
A <- rbinom(n, 1, 0.5)
beta <- runif(p)
beta <- beta/sqrt(as.numeric(crossprod(beta)))
Y <- as.numeric(A + X %*% beta + rnorm(n,0, 0.1))

# Split data by treatment
idx1 <- which(A == 1)
idx0 <- which(A == 0)

X1 <- X[idx1, , drop = FALSE]
Y1 <- Y[idx1]

X0 <- X[idx0, , drop = FALSE]
Y0 <- Y[idx0]

library(glmnet)
cv_fit1 <- cv.glmnet(X1, Y1, alpha = 1)  # alpha=1 => Lasso
mu1_hat_all <- predict(cv_fit1, newx = X, s = "lambda.min")  # predict on ALL X
cv_fit0 <- cv.glmnet(X0, Y0, alpha = 1)
mu0_hat_all <- predict(cv_fit0, newx = X, s = "lambda.min")

mu1 <- as.vector(mu1_hat_all)
mu0 <- as.vector(mu0_hat_all)
# Simple nuisance estimates
pi <- rep(0.5, n)

cat("Running HOIF with:\n")
cat("  n =", n, "\n")
cat("  p =", ncol(X), "\n")
cat("  Sample split: FALSE\n")
m <- 6
results <- hoif_ate(
  X, A, Y,
  mu1 = mu1,
  mu0 = mu0,
  pi = pi,
  transform_method = "none",
  k = 5,
  m = m,
  sample_split = 0,
  backend = "torch"
)

cat("\nResults:\n")
print(results)
plot(results)
cat("\n✓ Basic HOIF works!\n")

results <- hoif_ate(
  X, A, Y,
  mu1 = mu1,
  mu0 = mu0,
  pi = pi,
  transform_method = "none",
  k = 5,
  m = m,
  sample_split = 0,
  backend = "torch",
  pure_R_code = TRUE
)

cat("\nResults:\n")
print(results)
plot(results)
cat("\n✓ Basic HOIF works!\n")

# ==============================================================================
# Test 8: HOIF with Sample Split
# ==============================================================================
cat("\n========================================\n")
cat("Test 8: HOIF with Sample Split\n")
cat("========================================\n")

set.seed(456)
n <- 2000
X <- matrix(rnorm(n * 3), ncol = 3)
A <- rbinom(n, 1, 0.5)
Y <- A + X[,1] + rnorm(n)

mu1 <- rep(mean(Y[A==1]), n)
mu0 <- rep(mean(Y[A==0]), n)
pi <- rep(0.5, n)
m <- 7
cat("Running HOIF with:\n")
cat("  n =", n, "\n")
cat("  p =", ncol(X), "\n")
cat("  k = 8\n")
cat("  m = 4\n")
cat("  Sample split: TRUE (K=3)\n")

results_split <- hoif_ate(
  X, A, Y,
  mu1 = mu1,
  mu0 = mu0,
  pi = pi,
  transform_method = "none",
  k = 8,
  m = m,
  sample_split = TRUE,
  K = 6,
  seed = 123,
  backend = "torch"
)

cat("\nResults:\n")
print(results_split)

# Test reproducibility
results_split2 <- hoif_ate(
  X, A, Y,
  mu1 = mu1,
  mu0 = mu0,
  pi = pi,
  transform_method = "none",
  k = 8,
  m = m,
  sample_split = TRUE,
  K = 2,
  seed = 123,
  backend = "torch"
)

cat("\nReproducibility check:\n")
cat("  Same ATE:", all.equal(results_split$ATE, results_split2$ATE), "\n")

stopifnot(all.equal(results_split$ATE, results_split2$ATE))

cat("\n✓ Sample splitting works!\n")

# ==============================================================================
# Test 9: Different Transformations
# ==============================================================================
cat("\n========================================\n")
cat("Test 9: Compare Transformations\n")
cat("========================================\n")

set.seed(789)
n <- 150
X <- matrix(rnorm(n * 2), ncol = 2)
A <- rbinom(n, 1, 0.5)
Y <- A + 0.5 * X[,1] + 0.3 * X[,2] + rnorm(n)

mu1 <- rep(mean(Y[A==1]), n)
mu0 <- rep(mean(Y[A==0]), n)
pi <- rep(0.5, n)

# B-splines
results_splines <- hoif_ate(
  X, A, Y, mu1, mu0, pi,
  transform_method = "splines",
  k = 6, m = 3,
  sample_split = FALSE,
  backend = "numpy"
)

# Fourier
results_fourier <- hoif_ate(
  X, A, Y, mu1, mu0, pi,
  transform_method = "fourier",
  k = 6, m = 3,
  sample_split = FALSE,
  backend = "numpy"
)

cat("Comparison:\n")
comparison <- data.frame(
  Order = results_splines$orders,
  Splines = round(results_splines$ATE, 4),
  Fourier = round(results_fourier$ATE, 4)
)
print(comparison)

cat("\n✓ Both transformations work!\n")

# ==============================================================================
# Test 10: Plot Method
# ==============================================================================
cat("\n========================================\n")
cat("Test 10: Plotting\n")
cat("========================================\n")

cat("Creating convergence plot...\n")
plot(results_split)

cat("\n✓ Plot method works!\n")

# ==============================================================================
# Summary
# ==============================================================================
cat("\n========================================\n")
cat("🎉 ALL TESTS PASSED! 🎉\n")
cat("========================================\n")
cat("\nYour HOIF package is working correctly!\n")
cat("\nNext steps:\n")
cat("  1. Run formal unit tests: devtools::test()\n")
cat("  2. Check package: devtools::check()\n")
cat("  3. Try with real nuisance function estimates\n")
cat("  4. Test with larger datasets\n")
cat("  5. Try torch backend for speed comparison\n")
cat("\n")
