# R/utils_cdf.R

#' Install Custom CDF Package
#'
#' Installs a custom CDF package from mbni.org if not already installed.
#'
#' @param platform_name Cleaned platform name (e.g., "hugene10st").
#' @param species_code Two-letter species code (e.g., "hs", "mm").
#' @param cdf_version Version of the CDF (e.g., "25.0.0").
#' @param cdf_gene_map Gene mapping type (e.g., "ensg").
#' @return The name of the CDF package.
#' @importFrom utils install.packages
#' @keywords internal
install_custom_cdf <- function(platform_name, species_code, cdf_version, cdf_gene_map) {
  cdf_pkg_name <- paste0("pd", ".",
                         platform_name, ".",
                         species_code, ".",
                         cdf_gene_map)

  if (!requireNamespace(cdf_pkg_name, quietly = TRUE)) {
    message("Custom CDF package ", cdf_pkg_name, " not found. Attempting installation.")
    fp_pkg <- paste0(cdf_pkg_name, "_", cdf_version, ".tar.gz")
    url_pkg <- paste0("http://mbni.org/customcdf/",
                      cdf_version, "/", cdf_gene_map, ".download/",
                      basename(fp_pkg))

    # Check if the URL is reachable before attempting install.packages
    # This might require httr or similar, or a simple HEAD request check if possible.
    # For simplicity here, we proceed directly.
    message("Installing from: ", url_pkg)
    tryCatch({
      utils::install.packages(url_pkg, repos = NULL, type = "source")
      if (!requireNamespace(cdf_pkg_name, quietly = TRUE)) {
        stop("Failed to install and load custom CDF package: ", cdf_pkg_name)
      }
      message("Successfully installed ", cdf_pkg_name)
    }, error = function(e) {
      stop("Error installing custom CDF package ", cdf_pkg_name, " from ", url_pkg, ": ", e$message)
    })
  } else {
    message("Custom CDF package ", cdf_pkg_name, " is already installed.")
  }
  # Ensure the package is loaded for the current session if just installed
  # require() is generally discouraged in packages, but for dynamic package loading it might be okay.
  # Alternatively, the user might need to load it, or we rely on oligo to find it.
  # For oligo::read.celfiles, just having it installed is often enough if pkgname is correct.
  return(cdf_pkg_name)
}
