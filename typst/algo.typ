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
#let my-params = ("A", "Y", "X", "a_est", "b_est_1", "b_est_0", "order", "k", "Z_method", "is_split", "is_bootstrap", "bootstrap_seed", "bootstrap_number")

= main function
#algo(
  title: "HOIF",
  parameters: my-params
)[
 #smallcaps("Validate_parameters")[#my-params.join(", ")] \

  $n <-$ length of $A$ \
  $Z_k <-$ #smallcaps("Z_transform")$(k, X, "Z_method")$ \
  $epsilon_A <-$ #smallcaps("Epsilon_A")$(A, "a_est")$ \
  $epsilon_Y <-$ #smallcaps("Epsilon_Y")$(Y, A, "b_est_1", "b_est_0")$ \

  // if is_split is True:#i\
  //   $I_s <-$ sample indices for sample splitting\
  //   $Omega <-$ #smallcaps("Basis_Omega_estimation")($Z_k, A, I_s$) #d\
  // else:#i\
  //   $Omega <-$ #smallcaps("Basis_Omega_estimation")($Z_k, A$) #d\

  // $bold("IIFF") <-$ #smallcaps("IIFF_estimation")($Omega, Z_k, epsilon_A, epsilon_Y, "order", "is_split"$) \

  // if is_bootstrap is True:#i\
  //   Set random seed \
  //   Create weight matrix $W$ of size $n times n_"boot"$\
  //   for $i = 1$ to $n_"boot"$:#i\
  //     $bold("IIFF")_"boot"^{(i)} <-$ #smallcaps("IIFF_estimation") with weights $W_{*, i}$ #d\
  //   $V_"boot" <-$ #smallcaps("IIFF_var_estimation")($bold("IIFF")_"boot"$) #d\

  // return #smallcaps("Format_Results")($bold("IIFF"), V_"boot"$)
]

= Helper Functions



// #algo(
//   title: "Epsilon_A",
//   parameters: ("A", "a_est")
// )[
//   $epsilon_{A,1} <- 1 - A / a_est$\
//   $epsilon_{A,0} <- 1 - (1 - A) / (1 - a_est)$\
//   return $(epsilon_{A,1}, epsilon_{A,0})$
// ]

// #algo(
//   title: "Epsilon_Y",
//   parameters: ("Y", "A", "b_est_1", "b_est_0")
// )[
//   $epsilon_{Y,1} <- A dot (Y - b_est_1)$\
//   $epsilon_{Y,0} <- (1 - A) dot (Y - b_est_0)$\
//   return $(epsilon_{Y,1}, epsilon_{Y,0})$
// ]


// #algo(
//   title: "Basis_Omega_estimation",
//   parameters: ("Z_k", "A", "split_indices")
// )[
//   #comment[Estimate the projection basis or precision matrix]\
//   if split_indices is not null:#i\
//     Compute $Omega$ using only $Z_k$ and $A$ in split_indices #d\
//   else:#i\
//     Compute $Omega$ using full data #d\
//   return $Omega$
// ]

// #algo(
//   title: "IIFF_estimation",
//   parameters: ("Omega", "Z_k", "e_A", "e_Y", "order", "weights")
// )[
//   #comment[Placeholder for higher-order influence function calculation]\
//   Initialize $bold("IIFF")$ vector\
//   Apply $weights$ if provided\
//   Compute iterative expansion up to $order$\
//   return $bold("IIFF")$
// ]
