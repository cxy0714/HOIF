#import "@preview/clean-math-paper:0.2.4": *
#import "@preview/algo:0.3.6": algo, i, d, comment, code
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

$(X_i,A_i,Y_i)_(i=1)^(n)$ where $A_i in {0,1}, Y_i in bb(R), X_i in bb(R)^p$, we have the estimators of the nuisance functions $(hat(mu)(1,X_i), hat(mu)(0,X_i))$ and $pi(X_i)$, all function value is in $bb(R)$. That's the whole input data.

= Transformation of $X$

First we transform the covariates $X_i$ to some basis functions $Z_i in bb(R)^k$, B-splines or Fourier basis functions are recommended. 

here need a function 

input $X$ and transform method(fourier or splines) and the respective tuning parameters, k is also in the tuning parameter.

output: $(Z_i)_(i=1)^(n), Z_i in bb(R)^k$.

= Compute the residuals

We have 2 pair residuals, named by a label $a in {0,1}$.

when $a = 1$:

$ R_i^1 & = A_i (Y_i - hat(mu)(1,X_i)) \
  r_i^1 &  = 1 - A_i / hat(pi)(X_i)  $

When $a = 0$:
$ R_i^0 & = (1 - A_i) (Y_i - hat(mu)(0,X_i)) \
  r_i^0 & = 1 - (1 - A_i) / (1 - hat(pi)(X_i)) $

here may need a function 

input $(A_i,Y_i,hat(mu)(1,X_i),hat(mu)(0,X_i),hat(pi)(X_i))$ 

then output the 2 pair residuals $((R_i^1,r_i^1),(R_i^0,r_i^0))$.

= Compute the Inverse of weighted Gram matrix

we need to compute the inverse of weighted Gram matrix $G_a = 1/n sum_(i=1)^(n) s_i^a Z_i Z_i^T$ for $a in {0,1}$ where $s_i^a = A_i^a ( 1 - A_i)^(1-a)$, for compute the inverse matrix of $G_a$ named $Omega_a = G_a^(-1)$, we can use different method, direct inverse(using chol2inv() in R maybe fast) and the shrinkage method(nlshrink,corpcor package in R). 

here need a function

input $(Z_i, A_i)$ and the method(direct, nlshrink, corpcor), 

output $(Omega_1, Omega_0)$.

= Compute the basis matrix

compute the basis matrix for $a in {0,1}$:
$ B^a = Z Omega_a Z^T.  $
here $Z in bb(R)^(n times k)$, $Z_i$ is the i-th row of $Z$, 

here need a function

input $(Z_i)$ and $(Omega_1, Omega_0)$.
output $(B^1, B^0)$.

= Compute the HOIF estimators

finally we can compute the HOIF estimators for ATE:

input order $m$, residuals $((R_i^1,r_i^1),(R_i^0,r_i^0))$ and basis matrix $(B^1, B^0)$.

For order $j$ form $2$ to $m$, and for $a in {0,1}$, we compute the following statistics:
$ U_j^(a) = (-1)^j "ustat" ("tensors" = T_j^a, "expression"=  E_j^a, "backend" = "backend", "average" = 1). $
where $"ustat" $ is a pre-defined function to output a scalar, backend is input parameter valued in numpy or torch, and 
$ T_j^a &= "list"( R_i^a, underbrace(B^a, j-1"'s repeat"), r_i^a ) \
  E_j^a & = "list"(1, "list"(1,2), dots, "list"(j-1,j), j) $

Then for each order $l$ form $2$ to $m$ and  $a in (0,1)$ we compute the HOIF estimator:
$ "IIFF"_l^a &=  sum_(j=2)^(l) C_j^l U_j^(a) \
 "HOIF"_l^a & =  sum_(j=2)^(l) "IIFF"_j^a \
 "ATE"_l & = "HOIF"_l^1 - "HOIF"_l^0 $

where $C_j^l = binom(l-2,l-j)$. 

Then output 
$( "ATE"_l,"HOIF"_l^a,"IIFF"_l^a)$ for $l = 2, dots, m$ and $a in {0,1}$.


= The whole function - Sample spliting issue

The whole function is combining all above steps, input the whole data $(X_i,A_i,Y_i)_(i=1)^(n)$, the nuisance estimators $(hat(mu)(1,X_i), hat(mu)(0,X_i))$ and $hat(pi)(X_i)$, the transformation method and its tuning parameters, the inverse method of weighted Gram matrix, the order $m$ and the backend.

And we also need a switch to choose whether or not using sample spliting procedure and the number of splits $K$.

If not using sample spliting, above procedure is what we want, output  $( "ATE"_l,"HOIF"_l^a,"IIFF"_l^a)$ for $l = 2, dots, m$ and $a in {0,1}$ is over.

If using sample splitting, we need to split the whole data into $K$ parts first, say the indices are $(I_1, I_2, dots,I_k)$ conbines to whole indices ${1,2,dots,n}$, 

transformation of $X$ is done on the whole data first to get $(Z_i)_(i=1)^(n)$ and the residuals $((R_i^1,r_i^1),(R_i^0,r_i^0))$, but denote their $j$-th part with indices in $I_j$ as $(Z_i)_(i in I_j), ((R_i^1,r_i^1),(R_i^0,r_i^0))_(i in I_j)$.

Then

for j from 1 to K:

- use data  with indices not in $I_j$ i.e. ($(Z_i)_(i in.not I_j), (A_i)_(i in.not I_j))$ to compute the inverse of weighted Gram matrix $(Omega_(1,j), Omega_(0,j))$
- use data with indices in $I_j$ i.e. $(Z_i)_(i in I_j)$ and $(Omega_(1,j), Omega_(0,j))$ to compute the basis matrix $(B_(1,j), B_(0,j))$
- use data with indices in $I_j$ i.e. $((R_i^1,r_i^1),(R_i^0,r_i^0))_(i in I_j)$ and $(B_(1,j), B_(0,j))$ to compute the HOIF estimators $ ( "ATE"_(l,j), "HOIF"_(l,j)^a,"IIFF"_(l,j)^a) $ for $l = 2, dots, m$ and $a in {0,1}$.

end for

Finally average over $j$ from $1$ to $K$ to get the final estimators: 
$ "ATE"_l = 1/K sum_(j=1)^(K) "ATE"_(l,j) \
 "HOIF"_l^a = 1/K sum_(j=1)^(K) "HOIF"_(l,j)^a \
 "IIFF"_l^a = 1/K sum_(j=1)^(K) "IIFF"_(l,j)^a $
for $l = 2, dots, m$ and $a in {0,1}$.

This the whole function with sample spliting procedure.


#pagebreak()
= Implementation Notes for Developer (LLM Guidance)

#quote(block: true)[
  *Developer Instruction*: When implementing this algorithm in R (or Python), please pay strict attention to the following technical details to ensure the statistical validity of the HOIF estimators.
]

==  Basis Transformation Strategy
- *Input Scaling*: Always scale covariates $X$ to $[0, 1]$ before applying B-splines or Fourier transformations.
- *Additive Construction*: For $X in bb(R)^p$, generate bases for each dimension separately and concatenate them: $bold(Z) = [bold(1), phi(X_1), dots, phi(X_p)]$. Do not include interaction terms unless specified.
- *Intercept*: Ensure a constant term (column of 1s) is included in $bold(Z)$ to maintain the centering property of the Gram matrix.

==  Numerical Stability in Matrix Inversion
- *High-Dimensional Case*: If $k approx n$, the Gram matrix $bold(G)_a$ will be ill-conditioned. 
- *Efficiency*: In R, prefer `chol2inv(chol(Ga))` for speed, but wrap it in a `tryCatch` to handle non-positive-definite cases.

==  Sample Splitting (Cross-Fitting) Logic
This is the most critical part of the implementation. For each fold $k in {1, dots, K}$:
- *Training Set ($I_(-k)$)*: Used *only* to compute $bold(Omega)_{a, (-k)}$.
- *Estimation Set ($I_k$)*: Used to compute the local projection matrix $bold(B)_{a, k}$ and the $U$-statistics.
- *Index Alignment*: Ensure that the residuals $R_i^a$ and $r_i^a$ are indexed by $I_k$ and their order matches the rows/columns of $bold(B)_{a, k}$. 
  - $bold(B)_{a, k}$ should be of size $|I_k| times |I_k|$.
  - $R_("loc")$ and $r_("loc")$ should be vectors of length $|I_k|$.


==  Data Structure for Output
Return a nested list or a named list containing:
- `ATE`: A vector of length $m-1$ containing estimators from order $2$ to $m$.
- `IIFF_components`: The raw increment at each order to monitor convergence.
- `Convergence_Plot`: (Optional) A plot of $"ATE"_l$ vs. $l$ to verify if the estimator stabilizes as order increases.
