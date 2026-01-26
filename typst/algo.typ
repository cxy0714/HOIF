#import "@preview/clean-math-paper:0.2.4": *
#import "@preview/algo:0.3.6": algo, code, comment, d, i
#set par(first-line-indent: 1em)
#set page(margin: 1.75in)
#set par(leading: 0.55em, first-line-indent: 1.8em, justify: true)
#set text(font: "New Computer Modern")
#show par: set par(spacing: 0.55em)
#show heading: set block(above: 1.4em, below: 1em)

#let date = datetime.today().display("[month repr:long] [day], [year]")

// Modify some arguments, which can be overwritten in the template call
#page-args.insert("numbering", "1/1")
#text-args-title.insert("size", 2em)
#text-args-title.insert("fill", black)
#text-args-authors.insert("size", 12pt)

#show: template.with(
  title: "Algorithm of HOIF estimators for ATE",
  authors: (
    (name: "Xingyu Chen", affiliation-id: 1, orcid: "https://orcid.org/0009-0008-0823-4406"),
  ),
  affiliations: (
    (id: 1, name: "School of Mathematical Sciencei, Shanghai Jiao Tong University, China"),
  ),
  date: date,
  heading-color: rgb("#000000"),
  link-color: rgb("#008002"),
)



#outline()

= Input data

Let the observed data be denoted by $(X_i, A_i, Y_i)_{i=1}^{n}$, where $A_i in {0, 1}$ is the binary treatment indicator, $Y_i in bb(R)$ is the observed outcome, and $X_i in bb(R)^p$ represents the vector of covariates.

We assume the availability of pre-computed nuisance function estimators: the conditional mean outcomes $(hat(mu)(1, X_i), hat(mu)(0, X_i))$ and the propensity score $hat(pi)(X_i)$. All these estimators map to the real line $bb(R)$.



#let ustat_url = "https://github.com/your-repo/ustat"
#let ustat = link(ustat_url)[#math.upright("ustat")]

= The Integrated Algorithm (Main Interface)

The whole function combines all the following steps to get the HOIF estiamtors for ATE. It serves as the primary entry point, orchestrating the data transformation, residual calculation, and the choice of cross-fitting strategy.

*Function Inputs:*
- Full observed data $(X_i, A_i, Y_i)_(i=1)^n$.
- Nuisance function estimators $(hat(mu)(1, X_i), hat(mu)(0, X_i), hat(pi)(X_i))$.
- Transformation method and its tuning parameters.
- Inverse method of weighted Gram matrix.
- Maximum HOIF order $m$ and #ustat backend.
- **Switch**: A boolean flag for sample splitting and the number of splits $K$.

*Procedure:*

1. **Global Pre-processing**:
  - Transform the covariates $X_i$ on the whole data to obtain basis functions $(Z_i)_(i=1)^n$.
  - Compute the global residuals $((R_i^1, r_i^1), (R_i^0, r_i^0))_(i=1)^n$.

2. **Branching Logic**:
  - *If not using sample splitting*:
    Proceed with the entire dataset using the steps described in the following sections. Output $("ATE"_l, "HOIF"_l^a, "IIFF"_l^a)$ for $l = 2, dots, m$ and $a in {0,1}$.
  - *If using sample splitting (Cross-fitting)*:
    a. Split the indices ${1, dots, n}$ into $K$ disjoint parts $(I_1, I_2, dots, I_K)$.
    b. **For each fold** $j = 1, dots, K$:
    - **Training**: Use data with indices not in $I_j$ ($i in.not I_j$) to compute the inverse Gram matrices $(Omega_(1,j), Omega_(0,j))$.
    - **Estimation**: Use data with indices in $I_j$ ($i in I_j$) and the pre-computed $(Omega_(1,j), Omega_(0,j))$ to compute:
      - Local projection matrices $(B_(1,j), B_(0,j))$.
      - Local HOIF estimators $("ATE"_(l,j), "HOIF"_(l,j)^a, "IIFF"_(l,j)^a)$ for $l = 2, dots, m$.
    c. **Aggregation**: Average the results across all $K$ folds:
    $
      "ATE"_l = 1/K sum_(j=1)^(K) "ATE"_(l,j), quad "HOIF"_l^a = 1/K sum_(j=1)^(K) "HOIF"_(l,j)^a, quad "IIFF"_l^a = 1/K sum_(j=1)^(K) "IIFF"_(l,j)^a
    $

*Output:* Final averaged estimators $("ATE"_l, "HOIF"_l^a, "IIFF"_l^a)$ for $l = 2, dots, m$.

= Step-by-Step Technical Notes

== Transformation of Covariates

First, we transform the covariates $X_i$ into a set of basis functions $Z_i in bb(R)^k$. Common choices include B-splines or Fourier basis functions.

*Function Requirement:*
- *Input*: Covariate matrix $X$, the transformation method (e.g., Fourier or B-splines), and relevant tuning parameters (including the basis dimension $k$).
- *Output*: The transformed basis matrix $Z in bb(R)^(n times k)$, where each row $Z_i$ corresponds to the $i$-th observation.


== Compute the Residuals

Next, we compute two pairs of residuals, indexed by the treatment assignment $a in {0,1}$.

For the treatment group ($a = 1$):
$
  R_i^1 & = A_i (Y_i - hat(mu)(1,X_i)) \
  r_i^1 & = 1 - A_i / hat(pi)(X_i)
$

For the control group ($a = 0$):
$
  R_i^0 & = (1 - A_i) (Y_i - hat(mu)(0,X_i)) \
  r_i^0 & = 1 - (1 - A_i) / (1 - hat(pi)(X_i))
$

*Function Requirement:*
- *Input*: Observed data $(A, Y)$ and estimated nuisance functions $(hat(mu)(1,X), hat(mu)(0,X), hat(pi)(X))$.
- *Output*: Two residual pairs $((R_i^1, r_i^1))_(i=1)^n$ and $((R_i^0, r_i^0))_(i=1)^n$.


== Compute the Inverse of the Weighted Gram Matrix

Compute the inverse of the weighted Gram matrix $G_a$ for each $a in {0,1}$:
$ G_a = 1/n sum_(i=1)^(n) s_i^a Z_i Z_i^T, $
where $s_i^a = A_i^a (1 - A_i)^(1-a)$ serves as the indicator for the $a$-th group. Let $Omega_a = G_a^(-1)$ denote the corresponding inverse matrix.

Several estimation methods can be employed for $Omega_a$, such as direct inversion (e.g., using `chol2inv()` in R for efficiency) or shrinkage estimators (e.g., via the `nlshrink` or `corpcor` R packages) to improve numerical stability in high-dimensional settings.

*Function Requirement:*
- *Input*: Basis functions $Z$, treatment indicators $A$, and the inversion method (`direct`, `nlshrink`, or `corpcor`).
- *Output*: A pair of inverse matrices $(Omega_1, Omega_0)$.

== Compute the Projection Matrix

Compute the projection-like basis matrix $B^a$ for each $a in {0,1}$:
$ B^a = Z Omega_a Z^T, $
where $Z in bb(R)^(n times k)$ is the matrix of basis functions. Note that $B^a$ is an $n times n$ matrix.

*Function Requirement:*
- *Input*: Basis matrix $Z$ and the inverse matrices $(Omega_1, Omega_0)$.
- *Output*: Basis matrices $(B^1, B^0)$.


== Compute the HOIF Estimators

Finally, we compute the HOIF estimators for the ATE.

*Function Requirement:*
- *Input*:
  - Maximum order $m$.
  - Residual pairs $((R_i^a, r_i^a))_(a in {0,1})$.
  - Projection matrices $(B^1, B^0)$.
  - Backend for #ustat (`numpy` or `torch`).

- *Procedure*:
  For each order $j$ from $2$ to $m$, and for each treatment assignment $a in {0,1}$, calculate the $U$-statistics:
  $ "U"_j^(a) = (-1)^j #ustat ("tensors" = T_j^a, "expression" = E_j^a, "backend" = "backend", "average" = 1) $

  The tensors $T_j^a$ and index expressions $E_j^a$ are defined as:
  $
    T_j^a & = ( R^a, underbrace(B^a, j-1 " times"), r^a ) \
    E_j^a & = ( 1, (1,2), dots, (j-1, j), j )
  $

  Then, for each order $l in {2, dots, m}$, compute:
$
  "IIFF"_l^a & = sum_(j=2)^(l) binom(l-2, l-j) "U"_j^(a) \
  "HOIF"_l^a & = sum_(j=2)^(l) "IIFF"_j^a \
     "ATE"_l & = "HOIF"_l^1 - "HOIF"_l^0
$
- *Output*: $("ATE"_l, "HOIF"_l^a, "IIFF"_l^a)$ for $l = 2, dots, m$.

#pagebreak()
= Implementation Notes for Developer (LLM Guidance)

#quote(block: true)[
  *Developer Instruction*: When implementing this algorithm in R (or Python), please pay strict attention to the following technical details to ensure the statistical validity of the HOIF estimators.
]

== Basis Transformation Strategy
- *Input Scaling*: Always scale covariates $X$ to $[0, 1]$ before applying B-splines or Fourier transformations.
- *Additive Construction*: For $X in bb(R)^p$, generate bases for each dimension separately and concatenate them: $bold(Z) = [bold(1), phi(X_1), dots, phi(X_p)]$. Do not include interaction terms unless specified.
- *Intercept*: Ensure a constant term (column of 1s) is included in $bold(Z)$ to maintain the centering property of the Gram matrix.

== Numerical Stability in Matrix Inversion
- *High-Dimensional Case*: If $k approx n$, the Gram matrix $bold(G)_a$ will be ill-conditioned.
- *Efficiency*: In R, prefer `chol2inv(chol(Ga))` for speed, but wrap it in a `tryCatch` to handle non-positive-definite cases.

== Sample Splitting (Cross-Fitting) Logic
This is the most critical part of the implementation. For each fold $k in {1, dots, K}$:
- *Training Set ($I_(-k)$)*: Used *only* to compute $bold(Omega)_{a, (-k)}$.
- *Estimation Set ($I_k$)*: Used to compute the local projection matrix $bold(B)_{a, k}$ and the $U$-statistics.
- *Index Alignment*: Ensure that the residuals $R_i^a$ and $r_i^a$ are indexed by $I_k$ and their order matches the rows/columns of $bold(B)_{a, k}$.
  - $bold(B)_{a, k}$ should be of size $|I_k| times |I_k|$.
  - $R_("loc")$ and $r_("loc")$ should be vectors of length $|I_k|$.


== Data Structure for Output
Return a nested list or a named list containing:
- `ATE`: A vector of length $m-1$ containing estimators from order $2$ to $m$.
- `IIFF_components`: The raw increment at each order to monitor convergence.
- `Convergence_Plot`: (Optional) A plot of $"ATE"_l$ vs. $l$ to verify if the estimator stabilizes as order increases.
