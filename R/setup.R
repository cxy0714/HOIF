#' Environment Setup Functions for HOIF Package
#'
#' Functions for setting up and checking the Python environment required
#' for computing U-statistics.
#'
#' @author Xingyu Chen
#' @date 2026-01-23

#' Check if Python and u-stats are available
#'
#' @return Logical indicating availability
#' @keywords internal
check_python_env <- function() {
  # Check cached status first
  if (!is.null(hoif_env$python_available)) {
    return(hoif_env$python_available && hoif_env$ustats_available)
  }

  # Check Python
  py_ok <- tryCatch({
    reticulate::py_available(initialize = TRUE)
  }, error = function(e) FALSE)

  hoif_env$python_available <- py_ok

  if (!py_ok) {
    hoif_env$ustats_available <- FALSE
    return(FALSE)
  }

  # Check u-stats module
  ustats_ok <- tryCatch({
    reticulate::py_module_available("u_stats")
  }, error = function(e) FALSE)

  hoif_env$ustats_available <- ustats_ok

  return(py_ok && ustats_ok)
}


#' Setup Python environment for HOIF
#'
#' This function helps users set up the required Python environment
#' for using HOIF estimators. It will install Python (via Miniconda)
#' if needed and install the required u-stats package.
#'
#' @param method Character: installation method
#'   - "auto": automatically detect and use existing Python or install Miniconda
#'   - "virtualenv": create a virtual environment
#'   - "conda": use conda environment
#'   - "system": use system Python
#' @param envname Character: name of the virtual/conda environment (default: "r-hoif")
#' @param restart Logical: whether to restart R session after setup (default: FALSE)
#'
#' @return Invisible TRUE if successful
#'
#' @examples
#' \dontrun{
#' # First-time setup (recommended)
#' setup_hoif()
#'
#' # Use conda environment
#' setup_hoif(method = "conda", envname = "hoif-env")
#' }
#'
#' @export
setup_hoif <- function(method = c("auto", "virtualenv", "conda", "system"),
                       envname = "r-hoif",
                       restart = FALSE) {
  method <- match.arg(method)

  message("=== HOIF Python Environment Setup ===\n")

  # Step 1: Setup Python
  if (method == "auto") {
    if (reticulate::py_available(initialize = FALSE)) {
      message("✓ Found existing Python installation")
      py_config <- reticulate::py_config()
      message("  Python: ", py_config$python)
    } else {
      message("Installing Miniconda (this may take a few minutes)...")
      tryCatch({
        reticulate::install_miniconda()
        message("✓ Miniconda installed successfully")
      }, error = function(e) {
        stop("Failed to install Miniconda: ", e$message, call. = FALSE)
      })
    }
  } else if (method == "virtualenv") {
    message("Creating virtualenv: ", envname)
    reticulate::virtualenv_create(envname)
    reticulate::use_virtualenv(envname, required = TRUE)
    message("✓ Virtualenv created and activated")
  } else if (method == "conda") {
    message("Creating conda environment: ", envname)
    reticulate::conda_create(envname)
    reticulate::use_condaenv(envname, required = TRUE)
    message("✓ Conda environment created and activated")
  } else if (method == "system") {
    message("Using system Python")
    if (!reticulate::py_available(initialize = TRUE)) {
      stop("System Python not found. Please install Python first.", call. = FALSE)
    }
    message("✓ System Python detected")
  }

  # Step 2: Install required packages
  message("\nInstalling required Python packages...")

  required_pkgs <- c("u-stats", "numpy")

  tryCatch({
    reticulate::py_install(required_pkgs, pip = TRUE)
    message("✓ Required packages installed: ", paste(required_pkgs, collapse = ", "))
  }, error = function(e) {
    warning("Failed to install some packages: ", e$message, call. = FALSE)
    message("\nPlease install manually:")
    message("  pip install u-stats numpy")
    return(invisible(FALSE))
  })

  # Step 3: Optional PyTorch
  message("\nPyTorch can significantly speed up computations.")
  install_torch <- readline("Install PyTorch? (recommended) [y/N]: ")

  if (tolower(install_torch) %in% c("y", "yes")) {
    message("Installing PyTorch...")
    tryCatch({
      reticulate::py_install("torch", pip = TRUE)
      message("✓ PyTorch installed")
    }, error = function(e) {
      warning("Failed to install PyTorch: ", e$message, call. = FALSE)
      message("  You can install it later with: pip install torch")
    })
  }

  # Step 4: Verify installation
  message("\n=== Verifying Installation ===")

  # Reset cached status
  hoif_env$python_available <- NULL
  hoif_env$ustats_available <- NULL

  if (check_python_env()) {
    message("✓ Python environment is ready!")
    message("✓ u-stats package is available")

    # Test import
    tryCatch({
      ustats <- reticulate::import("u_stats")
      message("✓ Successfully imported u_stats")
      message("\n🎉 Setup complete! You can now use the HOIF package.")
    }, error = function(e) {
      warning("Could not import u_stats: ", e$message, call. = FALSE)
    })
  } else {
    warning("Setup completed but verification failed. You may need to restart R.",
            call. = FALSE)
    if (restart) {
      message("\nRestarting R session...")
      if (requireNamespace("rstudioapi", quietly = TRUE)) {
        rstudioapi::restartSession()
      } else {
        message("Please restart R manually (Ctrl+Shift+F10 in RStudio)")
      }
    } else {
      message("\nPlease restart your R session for changes to take effect.")
    }
  }

  invisible(TRUE)
}


#' Check HOIF environment status
#'
#' Checks whether the Python environment is properly configured for HOIF.
#'
#' @return Prints status message and invisibly returns TRUE/FALSE
#'
#' @examples
#' \dontrun{
#' check_hoif_setup()
#' }
#'
#' @export
check_hoif_setup <- function() {
  message("=== HOIF Environment Status ===\n")

  # Check Python
  py_ok <- tryCatch({
    reticulate::py_available(initialize = TRUE)
  }, error = function(e) FALSE)

  if (py_ok) {
    py_config <- reticulate::py_config()
    message("✓ Python: ", py_config$python)
    message("  Version: ", py_config$version)
  } else {
    message("✗ Python: Not found")
    message("\n  Run setup_hoif() to install Python")
    return(invisible(FALSE))
  }

  # Check u-stats
  ustats_ok <- tryCatch({
    reticulate::py_module_available("u_stats")
  }, error = function(e) FALSE)

  if (ustats_ok) {
    message("✓ u-stats: Available")
  } else {
    message("✗ u-stats: Not found")
    message("\n  Run setup_hoif() to install required packages")
    return(invisible(FALSE))
  }

  # Check torch (optional)
  torch_ok <- tryCatch({
    reticulate::py_module_available("torch")
  }, error = function(e) FALSE)

  if (torch_ok) {
    message("✓ PyTorch: Available (faster computation)")
  } else {
    message("○ PyTorch: Not found (optional, but recommended)")
    message("  Install with: reticulate::py_install('torch')")
  }

  # Check numpy
  numpy_ok <- tryCatch({
    reticulate::py_module_available("numpy")
  }, error = function(e) FALSE)

  if (numpy_ok) {
    message("✓ NumPy: Available")
  } else {
    message("✗ NumPy: Not found")
  }

  message("\n=================================")

  all_ok <- py_ok && ustats_ok && numpy_ok

  if (all_ok) {
    message("✓ Environment is ready!")
  } else {
    message("✗ Setup incomplete. Run setup_hoif() to fix issues.")
  }

  invisible(all_ok)
}
