# R/utils_geo.R

#' Generate GEO Prefix
#'
#' Creates a prefix string for constructing GEO FTP paths.
#'
#' @param geo_string A character string, typically a GSE accession number.
#' @return A character string representing the GEO prefix.
#' @keywords internal
#' @examples
#' # get_geo_prefix("GSE12345") # GSE12nnn
#' # get_geo_prefix("GSE123")   # GSEnnn
get_geo_prefix <- function(geo_string) {
  sapply(geo_string, function(bn) {
    if (nchar(bn) <= 6) {
      return(paste0(substr(bn, 1, 3),
                    paste(rep("n", 3),
                          collapse = "")))
    } else {
      return(paste0(substr(bn, 1, nchar(bn) - 3),
                    paste(rep("n", 3),
                          collapse = "")))
    }
  })
}

#' Download and Extract Raw GEO Data
#'
#' Downloads the _RAW.tar file for a given GSE and extracts it.
#'
#' @param gse_accession The GSE accession number (e.g., "GSE65656").
#' @param download_dir The directory where data should be downloaded and extracted.
#' @return The path to the directory containing the extracted raw files, or NULL if download fails.
#' @importFrom utils download.file untar
#' @keywords internal
download_and_extract_raw_geo <- function(gse_accession, download_dir) {
  raw_file_name <- paste0(gse_accession, "_RAW.tar")
  local_raw_file_path <- file.path(download_dir, raw_file_name)

  if (!file.exists(local_raw_file_path)) {
    remote_path_raw <- paste0("/geo/series/",
                              get_geo_prefix(gse_accession),
                              "/",
                              gse_accession,
                              "/suppl/",
                              raw_file_name)
    url_ftp <- paste0("ftp://ftp.ncbi.nlm.nih.gov", remote_path_raw)

    message("Downloading ", raw_file_name, " from ", url_ftp)
    tryCatch({
      utils::download.file(url = url_ftp, destfile = local_raw_file_path, mode = "wb", quiet = FALSE)
    }, error = function(e) {
      message("Failed to download ", raw_file_name, ": ", e$message)
      return(NULL)
    })
  }

  if (file.exists(local_raw_file_path)) {
    message("Extracting ", raw_file_name, " to ", download_dir)
    tryCatch({
      utils::untar(local_raw_file_path, exdir = download_dir)
      return(download_dir) # Or a more specific path if untar creates a subdir
    }, error = function(e) {
      message("Failed to extract ", raw_file_name, ": ", e$message)
      return(NULL)
    })
  } else {
    message("RAW file ", local_raw_file_path, " not found after download attempt.")
    return(NULL)
  }
}
