# zzz.R - Package loading hooks for HOIF

# .onAttach() - message on library(HOIF)
.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
    "HOIF: higher-order U-statistics are computed via the 'ustats' package.\n",
    "Its Python dependencies (u-stats, numpy, torch) are provisioned ",
    "automatically on first use; run ustats::check_ustats_setup() to verify ",
    "the environment, or see ?ustats::setup_ustats for manual setup.\n",
    "Alternatively, use hoif_ate(..., pure_R_code = TRUE) for a pure R ",
    "implementation (orders up to 6) that does not require Python."
  )
}
