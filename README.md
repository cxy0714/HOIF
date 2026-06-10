# HOIF

## Introduction

**HOIF** is an `R` package for the implementation of **Higher-Order Influence Function (HOIF)** estimators for the **Average Treatment Effect (ATE)**. The methodology is based on a series of foundational works by James M. Robins and his collaborators [1–4].

---

## The Core Computation

A core computational component of HOIF is the evaluation of **higher-order U-statistics**.

We have developed a general algorithm for computing U-statistics using the powerful Python functions [`numpy.einsum`](https://numpy.org/doc/stable/reference/generated/numpy.einsum.html) and [`torch.einsum`](https://pytorch.org/docs/stable/generated/torch.einsum.html).

We have built a Python package, [`u-stats`](https://pypi.org/project/u-stats/) (source: [U-Statistics-python](https://github.com/zrq1706/U-Statistics-python)), along with an R interface provided by [`ustats`](https://github.com/cxy0714/U-Statistics-R). Using [`ustats`](https://github.com/cxy0714/U-Statistics-R), we developed this `R` package for HOIF estimation of the ATE.

We also analyze computational complexity using graph-theoretic tools in a dedicated paper focused on the **exact computation** of higher-order U-statistics [5].

For HOIF estimators of the ATE, a key takeaway is:

- when the order $m \le 7$, the computational complexity is
  $O(n^3 + nk^2 + k^3 + n^2 k)$
- when the order $m > 7$, the computational complexity exceeds
  $O(n^4 + nk^2 + k^3 + n^2 k)$

Here, $n$ is the sample size and $k$ is the user-defined dimension of the transformed covariates $X$.
For more details on computing the higher-order U-statistics arising in HOIF, see Section 4.1 of [5].

The overall algorithmic workflow, mathematical formulas, and all parameters of `HOIF` are illustrated in [`inst/extdoc/HOIF.pdf`](inst/extdoc/HOIF.pdf) (after installation: `system.file("extdoc", "HOIF.pdf", package = "HOIF")`).

---

## Installation

```r
# From CRAN (once accepted):
# install.packages("HOIF")

# Development version from GitHub (also installs the dependency ustats):
# install.packages("devtools")
devtools::install_github("cxy0714/HOIF")
```

### Setting up the Python backend

By default, the higher-order U-statistics are computed in Python through the
[`ustats`](https://github.com/cxy0714/U-Statistics-R) package, which needs
**Python (>= 3.11)** with the Python packages `u-stats`, `numpy`, and
`torch`. There are three ways to get them — pick the one that fits you:

#### Option 1 — Do nothing (recommended)

With `reticulate` (>= 1.41), the Python dependencies are declared by `ustats`
and provisioned **automatically** the first time Python is needed. Just call
`hoif_ate()`:

```r
library(HOIF)
results <- hoif_ate(...)   # first call sets up Python automatically
```

The downloaded environment is cached and reused across sessions.

> **Note:** the first call downloads PyTorch. On Linux the default build
> bundles CUDA libraries (~2.5 GB); if you prefer a small CPU-only build
> or want full control, use Option 2 or 3.

#### Option 2 — One-shot managed setup: `ustats::setup_ustats()`

Creates a persistent environment and installs all dependencies. By default it
installs the **CPU-only** build of PyTorch (~200 MB instead of ~2.5 GB):

```r
ustats::setup_ustats()              # CPU-only PyTorch (default)
ustats::setup_ustats(gpu = TRUE)    # default PyPI PyTorch (CUDA on Linux)
```

#### Option 3 — Use your own Python/conda environment

If you already have an environment with a configured PyTorch (e.g. a specific
CUDA version), just add `u-stats`:

```bash
pip install u-stats
```

and point reticulate to that environment **before Python initializes**:

```r
library(HOIF)
reticulate::use_condaenv("your_env_name", required = TRUE)  # or use_virtualenv()
# Alternatively, set the RETICULATE_PYTHON environment variable
# (e.g. in .Rprofile or .Renviron) to the path of your python binary.
```

#### Verify the setup

Whichever option you chose:

```r
ustats::check_ustats_setup()
#> === ustats Environment Status ===
#> [OK] Python: /path/to/python
#> [OK] u_stats available
#> [OK] NumPy available
#> [OK] PyTorch available (version 2.x, CUDA available)
```

See `vignette("ustats", package = "ustats")` for a complete installation
guide and troubleshooting tips.

### No Python at all?

Pass `pure_R_code = TRUE` to `hoif_ate()`: the U-statistics are then computed
by a pure R implementation (orders up to 6) and **no Python runtime is
required**:

```r
results <- hoif_ate(X, A, Y, mu1 = mu1, mu0 = mu0, pi = pi,
                    m = 6, pure_R_code = TRUE)
```

---

## Example

A more comprehensive example can be found in [`test_manual/manual_test.R`](test_manual/manual_test.R).

```r
library(HOIF)

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


### Step 1: Estimate nuisance outcome regressions (here we omit the
### estimation of the propensity score by using the known randomization
### probability)

# - mu1(X) = E[Y | A = 1, X]
# - mu0(X) = E[Y | A = 0, X]

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
