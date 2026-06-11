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
  // authors: (
  //   (name: "Xingyu Chen", affiliation-id: 1, orcid: "https://orcid.org/0009-0008-0823-4406"),
  // ),
  // affiliations: (
  //   (id: 1, name: "School of Mathematical Sciencei, Shanghai Jiao Tong University, China"),
  // ),
  // date: date,
  heading-color: rgb("#000000"),
  link-color: rgb("#008002"),
)



#outline()


#let ustat_url = "https://github.com/cxy0714/U-Statistics-R"
#let ustat = link(ustat_url)[#math.upright("ustat")]

= Input data

Let the observed data be denoted by $(X_i, A_i, Y_i)_(i=1)^(n)$, where $A_i in {0, 1}$ is the binary treatment indicator, $Y_i in bb(R)$ is the observed outcome, and $X_i in bb(R)^p$ represents the vector of covariates.

We assume the availability of pre-computed nuisance function estimators: the conditional mean outcomes $(hat(mu)(1, X_i), hat(mu)(0, X_i))$ and the propensity score $hat(pi)(X_i)$. All these estimators map to the real line $bb(R)$. In standard applications, $hat(mu)(dot)$ and $hat(pi)(dot)$ should be estimated using an independent sample to avoid overfitting bias.


= Target formula

The core function `hoif_ate()` in this `R` package is to compute the below estimators.

Let $hat(psi)^("AIPW")_a$ denote the standard first-order doubly robust / double machine learning / AIPW estimator of $psi_a = bb(E)[Y(a)]$, and $hat("ATE")^("AIPW") = hat(psi)^("AIPW")_1 - hat(psi)^("AIPW")_0$. Conditional on the estimated nuisance functions, this estimator is biased.

$bb("BC")^("ATE")_m$ is the $m$-th order higher order influence function (HOIF) estimator of the *estimable bias* of the AIPW estimator of the average treatment effect (ATE) in causal inference, developed in a series of works by James M. Robins and his collaborators @Robins2008_HOIF @Robins2017_Minimax, @liu2017semiparametric, @LiuLi2023_sHOIF. It is *not* itself an estimator of the ATE; the sign convention is such that adding it to the AIPW estimate removes the estimable bias:
$
  hat("ATE")_m & = hat("ATE")^("AIPW") + bb("BC")^("ATE")_m \
$
(In the output of `hoif_ate()`, the components `HOIF1`, `HOIF0` and `ATE` correspond to $bb("HOIF")^1_m$, $bb("HOIF")^0_m$ and $bb("BC")^("ATE")_m$; the component name `ATE` is kept for backward compatibility -- it holds the ATE *bias correction*, not the ATE itself.)
$
  bb("BC")^("ATE")_m & = bb("HOIF")^1_m ( hat(Omega)^1) - bb("HOIF")^0_m ( hat(Omega)^0) \
  bb("HOIF")^a_m ( hat(Omega)^a) & = sum_(j=2)^m bb("IF")^a_j ( hat(Omega)^a) \
  bb("IF")^a_m ( hat(Omega)^a) & = (-1)^m ((n-m)!) /(n!) sum_((i_1, dots, i_m) in J_1^(times.o m) : i_1 eq.not i_2 eq.not dots eq.not i_m) r^a_(i_1) Z_(i_1)^(top) hat(Omega)^a product_(s = 2)^(m-1) {( Q^a_(i_s) -(hat(Omega)^a)^(-1) ) hat(Omega)^a} s^a_(i_m) Z_(i_m) R^a_(i_m) \
  & = sum_(j=2)^(m) binom(m-2, m-j) bb("U")^a_j ( hat(Omega)^a) \
  bb("U")^a_m ( hat(Omega)^a) & = (-1)^m ((n-m)!) /(n!) sum_((i_1, dots, i_m) in J_1^(times.o m) : i_1 eq.not i_2 eq.not dots eq.not i_m) r^a_(i_1) Z_(i_1)^(top) hat(Omega)^a product_(s = 2)^(m-1) { Q^a_(i_s) hat(Omega)^a} s^a_(i_m) Z_(i_m) R^a_(i_m) \
  & = (-1)^m ((n-m)!) /(n!) sum_((i_1, dots, i_m) in J_1^(times.o m) : i_1 eq.not i_2 eq.not dots eq.not i_m) r^a_(i_1) product_(s = 1)^(m-1) {Z_(i_s)^(top) hat(Omega)^a s^a_(i_(s+1)) Z_(i_(s+1)) } R^a_(i_m)
$

The notation above involves the following components:
$
       J_1, J_2 & subset.eq {1,2,dots,n}, \
              a & in {0,1}, \
        s^a_(i) & = A_i^a (1-A_i)^(1-a), \
        r^a_(i) & = 1 - s^a_(i) / ((hat(pi)(X_i))^a (1 - hat(pi)(X_i))^(1-a)), \
        R^a_(i) & = Y_i - hat(mu)(a,X_i), \
            Z_i & = serif("transform")(X_i) in bb(R)^k, \
        Q^a_(i) & = s^a_(i) Z_i Z_i^(top), \
   hat(Omega)^a & = (1/(|J_2|) sum_(i in J_2) s^a_(i) Z_i Z_i^(top))^(-1), \
$

Here, we designate the sample in $J_1$ as the _estimation sample_ and the sample in $J_2$ as the _training sample_ for $hat(Omega)^a$. The function $serif("transform")(X_i)$ denotes a transformation of the covariates, with the transformed variable $Z_i$ taking values in $bb(R)^k$.

*Sample Splitting (Cross-Fitting):* When employing $K$-fold sample splitting $(I_1, I_2, dots, I_K)$, for each fold $j in {1, dots, K}$, the sample in $I_j$ serves as the estimation sample (i.e., $I_j = J_1$), while samples not in $I_j$ form the training sample (i.e., $i in J_2 arrow.l.r.double i in.not I_j$). The training sample is used to compute $hat(Omega)^a$, which is then applied along with the estimation sample to calculate $bb("HOIF")^a_m$. The final estimator is obtained by averaging across all folds, yielding the empirical HOIF (eHOIF) estimator @liu2017semiparametric.

*No Sample Splitting:* When sample splitting is not used, we simply set $J_1 = J_2$, using the entire dataset for both training and estimation, which corresponds to the stable HOIF (sHOIF) estimator @LiuLi2023_sHOIF.


= The Integrated Algorithm (Main Interface)

The whole function combines all the following steps to get the HOIF bias-correction terms for the ATE. It serves as the primary entry point, orchestrating the data transformation, residual calculation, and the choice of cross-fitting strategy.

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
    Proceed with the entire dataset using the steps described in the following sections. Output $("BC"^("ATE")_l, "HOIF"_l^a, "IIFF"_l^a)$ for $l = 2, dots, m$ and $a in {0,1}$ (returned by the package as `ATE`, `HOIF1`/`HOIF0`, `IIFF1`/`IIFF0`).
  - *If using sample splitting (Cross-fitting)*:
    a. Split the indices ${1, dots, n}$ into $K$ disjoint parts $(I_1, I_2, dots, I_K)$.
    b. **For each fold** $j = 1, dots, K$:
    - **Training**: Use data with indices not in $I_j$ ($i in.not I_j$) to compute the inverse Gram matrices $(Omega_(1,j), Omega_(0,j))$.
    - **Estimation**: Use data with indices in $I_j$ ($i in I_j$) and the pre-computed $(Omega_(1,j), Omega_(0,j))$ to compute:
      - Local projection matrices $(B_(1,j), B_(0,j))$.
      - Local HOIF bias-correction terms $("BC"^("ATE")_(l,j), "HOIF"_(l,j)^a, "IIFF"_(l,j)^a)$ for $l = 2, dots, m$.
    c. **Aggregation**: Average the results across all $K$ folds:
    $
      "BC"^("ATE")_l = 1/K sum_(j=1)^(K) "BC"^("ATE")_(l,j), quad "HOIF"_l^a = 1/K sum_(j=1)^(K) "HOIF"_(l,j)^a, quad "IIFF"_l^a = 1/K sum_(j=1)^(K) "IIFF"_(l,j)^a
    $

*Output:* Final averaged bias-correction terms $("BC"^("ATE")_l, "HOIF"_l^a, "IIFF"_l^a)$ for $l = 2, dots, m$.

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
  R_i^1 & = Y_i - hat(mu)(1,X_i) \
  r_i^1 & = 1 - A_i / hat(pi)(X_i)
$

For the control group ($a = 0$):
$
  R_i^0 & = Y_i - hat(mu)(0,X_i) \
  r_i^0 & = 1 - (1 - A_i) / (1 - hat(pi)(X_i))
$

(The treatment indicator $s^a_i$ is *not* folded into $R^a_i$ here; it enters through the projection matrix $B^a$ below, matching the implementation.)

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
$ B^a = Z Omega_a (s^a dot.circle Z)^T, quad "i.e." quad B^a_(i j) = Z_i^top Omega_a s^a_j Z_j, $
where $Z in bb(R)^(n times k)$ is the matrix of basis functions and $s^a dot.circle Z$ multiplies row $j$ of $Z$ by $s^a_j$ (in code: `B1 = Z %*% Omega1 %*% t(Z * A)`). Note that $B^a$ is an $n times n$ matrix carrying the treatment indicator on its second index.

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
    T_j^a & = ( r^a, underbrace(B^a, j-1 " times"), R^a ) \
    E_j^a & = ( 1, (1,2), dots, (j-1, j), j )
  $
  so that the first tensor $r^a$ is indexed by $i_1$ and the last tensor $R^a$ by $i_j$, matching the target formula.

  Then, for each order $l in {2, dots, m}$, compute:
$
  "IIFF"_l^a & = sum_(j=2)^(l) binom(l-2, l-j) "U"_j^(a) \
  "HOIF"_l^a & = sum_(j=2)^(l) "IIFF"_j^a \
     "BC"^("ATE")_l & = "HOIF"_l^1 - "HOIF"_l^0
$
- *Output*: $("BC"^("ATE")_l, "HOIF"_l^a, "IIFF"_l^a)$ for $l = 2, dots, m$, returned by the package as `ATE`, `HOIF1`/`HOIF0` and `IIFF1`/`IIFF0`.




#bibliography("Master.bib")
