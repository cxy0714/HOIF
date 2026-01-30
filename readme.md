# HOIF

## Introduction

**HOIF** is an `R` package for the implementation of **Higher-Order Influence Function (HOIF)** estimators for the **Average Treatment Effect (ATE)**. The methodology is based on a series of foundational works by James M. Robins and his collaborators [1–4].

---

## The Core Computation

A core computational component of HOIF is the evaluation of **higher-order U-statistics**.  

We have developed a general algorithm for computing  U-statistics using the powerful Python functions [`numpy.einsum`](https://numpy.org/doc/stable/reference/generated/numpy.einsum.html) and [`torch.einsum`](https://pytorch.org/docs/stable/generated/torch.einsum.html).

We have built a Python package, [`ustats-python`](https://github.com/zrq1706/U-Statistics-python/tree/main), along with an R interface provided by [`ustats-R`](https://github.com/cxy0714/U-Statistics-R). Using [`ustats-R`](https://github.com/cxy0714/U-Statistics-R), we developed this `R` package for HOIF estimation of the ATE.

We also analyze computational complexity using graph-theoretic tools in a dedicated paper focused on the **exact computation** of higher-order U-statistics [5].

For HOIF estimators of the ATE, a key takeaway is:

- when the order $m \le 7$, the computational complexity is  
  $O(n^3 + nk^2 + k^3 + n^2 k)$  
- when the order $m > 7$, the computational complexity exceeds  
  $O(n^4 + nk^2 + k^3 + n^2 k)$

Here, $n$ is the sample size and $k$ is the user-defined dimension of the transformed covariates $X$.  
For more details on computing the higher-order U-statistics arising in HOIF, see Section 4.1 of [5].

The overall algorithmic workflow, mathematical formulas, and all parameters of `HOIF` are illustrated in [`inst/latex/HOIF.pdf`](inst/latex/HOIF.pdf).


---

## Installation

Install the development version from GitHub (this will also install the dependency [`ustats-R`](https://github.com/cxy0714/U-Statistics-R)):

```r
devtools::install_github("cxy0714/HOIF")
```

Load the required libraries and set up the `reticulate` environment using `setup_ustats()` from [`ustats-R`](https://github.com/cxy0714/U-Statistics-R):


```r
library(ustats)
library(HOIF)

library(reticulate)
py_config()

setup_ustats()
check_ustats_setup()
```

---

## Example

A more comprehensive example can be found in [`test/manual_test.R`](test_manual/manual_test.R).

```r
set.seed(123)

n <- 1000   # Sample size
p <- 50     # Number of covariates

# Covariates
X <- matrix(rnorm(n * p), ncol = p)

# Binary treatment assignment (randomized)
A <- rbinom(n, 1, 0.5)

# True coefficient vector (normalized)
beta <- runif(p)
beta <- beta / sqrt(as.numeric(crossprod(beta)))

# Outcome: partially linear model
Y <- as.numeric(A + X %*% beta + rnorm(n, 0, 0.1))


### Step 1: Estimate nuisance outcome regressions here omit the estimate of propensity score function by set the true propensity score

- μ₁(X) = E[Y | A = 1, X]  
- μ₀(X) = E[Y | A = 0, X]

idx1 <- which(A == 1)
idx0 <- which(A == 0)

X1 <- X[idx1, , drop = FALSE]
Y1 <- Y[idx1]

X0 <- X[idx0, , drop = FALSE]
Y0 <- Y[idx0]

library(glmnet)

cv_fit1 <- cv.glmnet(X1, Y1, alpha = 1)
mu1_hat_all <- predict(cv_fit1, newx = X, s = "lambda.min")

cv_fit0 <- cv.glmnet(X0, Y0, alpha = 1)
mu0_hat_all <- predict(cv_fit0, newx = X, s = "lambda.min")

mu1 <- as.vector(mu1_hat_all)
mu0 <- as.vector(mu0_hat_all)

# Known propensity score
pi <- rep(0.5, n)


### Step 2: Run HOIF with sample splitting (eHOIF)


m <- 7
n_folds <- 2

results_split <- hoif_ate(
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
print(results_split)
plot(results_split)

### Step 3: Run HOIF without sample splitting (sHOIF)


m <- 7

results_no_split <- hoif_ate(
  X, A, Y,
  mu1 = mu1,
  mu0 = mu0,
  pi = pi,
  transform_method = "none",
  m = m,
  sample_split = FALSE,
  backend = "torch"
)

cat("\nResults (without sample splitting):\n")
print(results_no_split)
plot(results_no_split)
```

---

## References

<sup>1</sup> Robins, J., Li, L., Tchetgen Tchetgen, E., & van der Vaart, A. (2008). Higher order influence functions and minimax estimation of nonlinear functionals. In Probability and Statistics: Essays in Honor of David A. Freedman, 335–421.

<sup>2</sup> Robins, J. M., Li, L., Mukherjee, R., Tchetgen Tchetgen, E., & van der Vaart, A. (2017). Minimax estimation of a functional on a structured high-dimensional model. The Annals of Statistics, 45(5), 1951–1987.

<sup>3</sup> Liu, L., Mukherjee, R., Newey, W. K., & Robins, J. M. (2017). Semiparametric efficient empirical higher order influence function estimators. arXiv:1705.07577.

<sup>4</sup> Liu, L., & Li, C. (2023). New √n-consistent, numerically stable higher-order influence function estimators. arXiv:2302.08097.

<sup>5</sup> Chen, X., Zhang, R., & Liu, L. (2025). On computing and the complexity of computing higher-order U-statistics, exactly. arXiv:2508.12627.
