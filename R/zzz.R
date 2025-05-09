# R/zzz.R

# These could also be arguments to the main function with these as defaults.
# Storing them here makes them accessible within the package but not exported.
.ExpressoGEOEnv <- new.env(parent = emptyenv())

.ExpressoGEOEnv$DEFAULT_CDF_VERSION <- "25.0.0"
.ExpressoGEOEnv$DEFAULT_CDF_GENE_MAP <- "ensg" # ensembl 102

.ExpressoGEOEnv$SPECIES_MAP <- c(
  "Mus musculus" = "mm",
  "Homo sapiens" = "hs",
  "Rattus norvegicus" = "rn"
)

.onLoad <- function(libname, pkgname) {
  required_bioc_pkgs <- c(
    "GEOquery", "Biobase", "affyio", "affy",
    "oligo", "AnnotationHub", "ensembldb"
  )

  missing_pkgs <- character(0)
  for (pkg in required_bioc_pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      missing_pkgs <- c(missing_pkgs, pkg)
    }
  }

  if (length(missing_pkgs) > 0) {
    missing_msg <- paste(
      "The following essential Bioconductor packages are not installed: ",
      paste(missing_pkgs, collapse = ", "),
      ". Please install them by running: \n\n",
      "if (!requireNamespace(\"BiocManager\", quietly = TRUE))\n",
      "    install.packages(\"BiocManager\")\n",
      "BiocManager::install(c('", paste(missing_pkgs, collapse = "', '"), "'))\n\n",
      "Then, try loading ", pkgname, " again.",
      sep = ""
    )
    # packageStartupMessage is preferred over stop() in .onLoad for non-critical load failures
    # that just prevent functionality. If the package is unusable, stop() is okay.
    packageStartupMessage(warning(missing_msg, call. = FALSE))
    # Or if absolutely critical:
    # stop(missing_msg, call. = FALSE)
  }
}
# .onUnload <- function(libpath) {
#   # Code to run when package is unloaded
# }
