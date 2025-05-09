# R/process_gse.R

#' Process Affymetrix Microarray Data from GEO
#'
#' Downloads, processes, and annotates Affymetrix microarray data for a given
#' Gene Expression Omnibus (GEO) Series (GSE) accession number.
#' It handles downloading raw CEL files, installing custom CDF packages if needed,
#' performing RMA normalization, and annotating genes using Ensembl via AnnotationHub.
#'
#' @param gse_accession A character string specifying the GSE accession number (e.g., "GSE65656").
#' @param output_dir A character string specifying the directory where downloaded files
#'        and the final RDS output will be saved. Defaults to a temporary directory.
#'        If NULL, RDS is not saved, and ExpressionSet is returned.
#' @param cdf_version A character string for the custom CDF version. Defaults to "25.0.0".
#' @param cdf_gene_map A character string for the gene mapping type in custom CDF. Defaults to "ensg".
#' @param force_reprocess Logical, if TRUE, will reprocess even if an output RDS file exists.
#'        Defaults to FALSE.
#'
#' @return An \code{ExpressionSet} object containing the normalized and annotated data.
#'         If \code{output_dir} is specified, the \code{ExpressionSet} is also saved as an RDS file
#'         in that directory. Returns NULL if processing fails for a platform.
#'
#' @details
#' This function automates several steps:
#' \enumerate{
#'   \item Validates the GSE accession format.
#'   \item Fetches metadata for the GSE using \code{GEOquery}.
#'   \item Downloads and extracts the _RAW.tar archive containing CEL files.
#'   \item For each platform within the GSE:
#'     \itemize{
#'       \item Determines the species and microarray platform type.
#'       \item If species and platform are unique for the samples on that platform:
#'         \itemize{
#'           \item Constructs and installs a custom CDF package from \url{http://mbni.org/customcdf/} if not already present.
#'           \item Reads CEL files using \code{oligo::read.celfiles} with the custom CDF.
#'           \item Performs RMA normalization using \code{oligo::rma}.
#'           \item Annotates probe IDs to Ensembl gene IDs and gene names using \code{AnnotationHub}.
#'           \item Creates and returns an \code{ExpressionSet} object.
#'           \item Saves the \code{ExpressionSet} to an RDS file if \code{output_dir} is provided.
#'         }
#'       \item Skips processing for platforms with mixed species or microarray types.
#'     }
#' }
#'
#' @examples
#' \dontrun{
#' # Ensure you have a working internet connection and write permissions
#' # to a directory if not using tempdir().
#' # This example uses a known small GSE for testing.
#' # Replace "GSE19804" with a GSE of interest if it uses a supported platform.
#' # Note: GSE19804 uses hgu133plus2, which might not have a direct mbni.org ensg CDF.
#' # You might need to find a GSE that uses platforms like MoGene or HuGene arrays
#' # for the custom CDF part to work as intended by the original script.
#'
#' # Example with a hypothetical GSE known to work with custom CDFs:
#' # gse_id <- "GSE_EXAMPLE_CUSTOM_CDF" # Replace with a real one
#' # eset <- process_gse(gse_id, output_dir = "./processed_data")
#' # if (!is.null(eset)) {
#' #   print(eset)
#' # }
#'
#' # To run with a temporary directory:
#' # eset_temp <- process_gse("GSE_EXAMPLE_CUSTOM_CDF") # Replace
#' }
#'
#' @seealso \code{\link[GEOquery]{getGEO}}, \code{\link[oligo]{read.celfiles}}, \code{\link[oligo]{rma}}, \code{\link[AnnotationHub]{AnnotationHub}}, \code{\link[Biobase]{ExpressionSet}}
#'
#' @export
#' @importFrom GEOquery getGEO
#' @importFrom Biobase pData exprs ExpressionSet AnnotatedDataFrame featureNames featureNames<- sampleNames sampleNames<-
#' @importFrom affyio read.celfile.header
#' @importFrom affy cleancdfname
#' @importFrom oligo read.celfiles rma
#'
process_gse <- function(gse_accession, output_dir = tempdir(),
                        cdf_version = "25.0.0", cdf_gene_map = "ensg",
                        force_reprocess = FALSE) {

  # --- Argument Validation ---
  if (!is.character(gse_accession) || length(gse_accession) != 1 || !grepl("^GSE\\d+$", gse_accession)) {
    stop("Invalid GSE accession number. Must be a string like 'GSE12345'.")
  }
  if (!is.null(output_dir) && (!is.character(output_dir) || length(output_dir) != 1)) {
    stop("'output_dir' must be a single character string or NULL.")
  }
  # Add more validation for cdf_version, cdf_gene_map if necessary

  # --- Setup Constants (could be from zzz.R or package environment) ---
  species_name_to_code_map <- c(
    "Mus musculus" = "mm",
    "Homo sapiens" = "hs",
    "Rattus norvegicus" = "rn"
  )

  # --- Directory Setup ---
  if (!is.null(output_dir)) {
    gse_data_dir <- file.path(output_dir, gse_accession)
    if (!dir.exists(gse_data_dir)) {
      dir.create(gse_data_dir, recursive = TRUE)
      message("Created directory: ", gse_data_dir)
    }
  } else {
    # If output_dir is NULL, use a temporary directory for downloads that gets cleaned up
    gse_data_dir <- tempfile(pattern = paste0(gse_accession, "_data_"))
    dir.create(gse_data_dir, recursive = TRUE)
    on.exit(unlink(gse_data_dir, recursive = TRUE, force = TRUE), add = TRUE)
    message("Using temporary directory for downloads: ", gse_data_dir)
  }


  message("Processing GSE: ", gse_accession)

  # --- Get GEO Metadata ---
  geo_metadata_list <- tryCatch({
    GEOquery::getGEO(
      gse_accession,
      destdir = gse_data_dir, # Downloads GPL, series matrix etc. here
      GSEMatrix = TRUE,
      AnnotGPL = FALSE, # We handle annotation later
      getGPL = FALSE    # Usually not needed if not AnnotGPL
    )
  }, error = function(e) {
    stop("Failed to retrieve GEO metadata for ", gse_accession, ": ", e$message)
  })

  # --- Download and Extract Raw Data ---
  # This will download to gse_data_dir/GSExxxx_RAW.tar and extract there
  raw_data_path <- download_and_extract_raw_geo(gse_accession, gse_data_dir)
  if (is.null(raw_data_path)) {
    warning("Failed to download or extract RAW data for ", gse_accession, ". Skipping.")
    return(NULL)
  }

  # --- Process each platform/series in the GSE ---
  # A GSE can contain multiple sub-series (GPLs)
  all_esets <- list() # To store ExpressionSet for each successfully processed platform

  for (i_plt in seq_along(geo_metadata_list)) {
    gset_platform <- geo_metadata_list[[i_plt]]
    platform_id <- Biobase::annotation(gset_platform) # GPL id
    message("\nProcessing platform: ", platform_id, " (", i_plt, " of ", length(geo_metadata_list), ")")

    pheno_data_df <- Biobase::pData(gset_platform)

    # Check for supplementary file paths (CEL files)
    if (!"supplementary_file" %in% colnames(pheno_data_df) ||
        all(sapply(pheno_data_df$supplementary_file, function(x) all(is.na(x) | x == "")))) {
      message("No supplementary files (CEL files) listed for platform ", platform_id, ". Skipping.")
      next
    }

    # Ensure supplementary files are full paths
    # getGEO sometimes returns relative paths based on its destdir
    # The raw data has been extracted into gse_data_dir
    # Basenames of supplementary files should now exist in gse_data_dir
    pheno_data_df$supplementary_file_full_path <- file.path(
      gse_data_dir,
      basename(as.character(pheno_data_df$supplementary_file))
    )

    # Check if CEL files actually exist
    cel_files_exist <- file.exists(pheno_data_df$supplementary_file_full_path)
    if (!all(cel_files_exist)) {
      message("Missing some CEL files for platform ", platform_id, ". Expected in: ", gse_data_dir)
      message("Missing: ", paste(basename(pheno_data_df$supplementary_file_full_path[!cel_files_exist]), collapse=", "))
      # Attempt to find any .CEL.gz or .cel.gz files in gse_data_dir
      all_extracted_cels <- list.files(gse_data_dir, pattern = "\\.cel(\\.gz)?$", ignore.case = TRUE, full.names = TRUE, recursive = TRUE)
      if(length(all_extracted_cels) < nrow(pheno_data_df)){
        message("Found only ", length(all_extracted_cels), " CEL files in ", gse_data_dir, " but expected ", nrow(pheno_data_df), ". Skipping platform ", platform_id)
        next
      } else {
        # This part is tricky if filenames don't match supplementary_file column.
        # For now, we rely on supplementary_file column being correct.
        # If your workflow often has mismatched names, this needs more robust logic.
        message("Strict check on supplementary_file names failed. If issues persist, verify CEL file names and paths.")
      }
    }


    # --- Determine Species ---
    # organism_ch1 can be a list if samples have different organisms.
    # We require a single organism for this processing pipeline.
    organism_values <- unique(unlist(lapply(pheno_data_df$organism_ch1, function(x) strsplit(as.character(x), ";")[[1]])))
    organism_values <- trimws(organism_values) # Trim whitespace

    if (length(organism_values) > 1) {
      message("Multiple organisms found for platform ", platform_id, ": ", paste(organism_values, collapse=", "), ". Skipping.")
      next
    }
    gse_species_name <- organism_values[1]
    gse_species_code <- species_name_to_code_map[gse_species_name]

    if (is.na(gse_species_code)) {
      message("Unsupported organism: ", gse_species_name, " for platform ", platform_id, ". Skipping.")
      next
    }
    message("Organism: ", gse_species_name, " (Code: ", gse_species_code, ")")

    # --- Determine Microarray Platform Type from CEL headers ---
    platform_names_from_cel <- sapply(pheno_data_df$supplementary_file_full_path, function(suppl_file_path) {
      if (!file.exists(suppl_file_path)) {
        warning("CEL file not found: ", suppl_file_path)
        return(NA_character_)
      }
      tryCatch({
        header_info <- affyio::read.celfile.header(suppl_file_path, info = "full")
        cleaned_cdf_name <- affy::cleancdfname(header_info$cdfName, addcdf = FALSE)
        # Remove potential version suffixes like "v1", "v2" for broader matching
        cleaned_cdf_name <- gsub("v\\d+$", "", cleaned_cdf_name)
        return(cleaned_cdf_name)
      }, error = function(e) {
        warning("Could not read CEL header for: ", basename(suppl_file_path), " - ", e$message)
        return(NA_character_)
      })
    })

    platform_names_from_cel <- platform_names_from_cel[!is.na(platform_names_from_cel)]
    unique_platform_types <- unique(platform_names_from_cel)

    if (length(unique_platform_types) == 0) {
      message("Could not determine platform type from any CEL files for platform ", platform_id, ". Skipping.")
      next
    }
    if (length(unique_platform_types) > 1) {
      message("Multiple microarray platform types found within samples for GPL ", platform_id, ": ",
              paste(unique_platform_types, collapse=", "), ". Skipping this set.")
      next
    }
    gse_platform_type <- unique_platform_types[1]
    message("Microarray Platform Type: ", gse_platform_type)

    # --- Define Output File Name ---
    output_rds_filename <- NULL
    if (!is.null(output_dir)) {
      output_rds_filename <- file.path(
        output_dir, # Main output_dir provided by user
        paste0(gse_species_code, "_", gse_platform_type, "_", gse_accession, "_", platform_id, ".RDS")
      )
      if (file.exists(output_rds_filename) && !force_reprocess) {
        message("Processed data already exists: ", output_rds_filename, ". Skipping. Use force_reprocess=TRUE to override.")
        # Optionally load and return existing RDS
        # all_esets[[platform_id]] <- readRDS(output_rds_filename)
        next
      }
    }


    # --- Install Custom CDF if necessary ---
    message("Attempting to set up CDF package for platform: ", gse_platform_type, ", species: ", gse_species_code)
    cdf_package_name <- tryCatch({
      install_custom_cdf(gse_platform_type, gse_species_code, cdf_version, cdf_gene_map)
    }, error = function(e) {
      message("Failed to install or verify custom CDF for ", gse_platform_type, " / ", gse_species_code, ": ", e$message)
      message("This platform might not be supported by the mbni.org custom CDFs or there was an installation issue.")
      message("You may need to install a suitable pd.<platform>.<species>.<map> package manually or use a different annotation method.")
      return(NULL)
    })

    if (is.null(cdf_package_name)) {
      message("Skipping processing for platform ", platform_id, " due to CDF issues.")
      next
    }

    # --- Read CEL files and Normalize ---
    message("Reading CEL files using oligo and CDF: ", cdf_package_name)
    raw_oligo_data <- tryCatch({
      oligo::read.celfiles(
        filenames = pheno_data_df$supplementary_file_full_path,
        sampleNames = pheno_data_df$geo_accession, # Use GSM IDs as sample names
        phenoData = new("AnnotatedDataFrame", data = pheno_data_df), # Pass full pData
        pkgname = cdf_package_name
      )
    }, error = function(e) {
      message("Error reading CEL files with oligo for platform ", platform_id, ": ", e$message)
      return(NULL)
    })

    if (is.null(raw_oligo_data)) {
      message("Skipping normalization for platform ", platform_id, ".")
      next
    }

    message("Performing RMA normalization...")
    normalized_eset <- tryCatch({
      oligo::rma(raw_oligo_data)
    }, error = function(e) {
      message("Error during RMA normalization for platform ", platform_id, ": ", e$message)
      return(NULL)
    })

    if (is.null(normalized_eset)) {
      message("Skipping further processing for platform ", platform_id, ".")
      next
    }

    # Clean rownames (remove "_at" common in some Affy arrays, though custom CDFs might use Ensembl IDs directly)
    Biobase::featureNames(normalized_eset) <- gsub("_at$", "", Biobase::featureNames(normalized_eset))
    current_gene_ids <- Biobase::featureNames(normalized_eset)

    # --- Annotate Features ---
    message("Annotating features using AnnotationHub for species: ", gse_species_code)
    feature_annotation_df <- tryCatch({
      fetch_gene_annotations(gse_species_code, current_gene_ids)
    }, error = function(e) {
      message("Error fetching gene annotations for platform ", platform_id, ": ", e$message)
      message("ExpressionSet will be created with available feature IDs but limited gene names.")
      # Create a basic feature data frame if annotation fails
      data.frame(gene_id = current_gene_ids, row.names = current_gene_ids)
    })

    # Ensure the featureData has rownames matching the ExpressionSet featureNames
    # fetch_gene_annotations should already set rownames to gene_id
    # If feature_annotation_df row order or content doesn't perfectly match normalized_eset,
    # this step will align them.
    # Create an empty data frame with correct rownames if annotation failed completely
    if (nrow(feature_annotation_df) == 0 || !identical(sort(rownames(feature_annotation_df)), sort(current_gene_ids))) {
      warning("Mismatch or failure in feature annotation. Using basic feature IDs.")
      final_feature_data_df <- data.frame(gene_id = current_gene_ids, row.names = current_gene_ids)
    } else {
      # Reorder feature_annotation_df to match current_gene_ids (rownames of exprs matrix)
      final_feature_data_df <- feature_annotation_df[current_gene_ids, , drop = FALSE]
      rownames(final_feature_data_df) <- current_gene_ids # Ensure rownames are set
    }

    # --- Create Final ExpressionSet ---
    # PhenoData comes from the normalized_eset (which got it from raw_oligo_data, which got it from gset_platform)
    # Or, ensure it's directly from gset_platform, re-aligned if necessary
    # The sampleNames of normalized_eset should match row.names of pData(gset_platform)
    # If oligo changed sample names, this needs care. `oligo::read.celfiles` used geo_accession.
    # So pData(normalized_eset) should be correct.

    final_eset <- Biobase::ExpressionSet(
      assayData = Biobase::exprs(normalized_eset),
      phenoData = Biobase::phenoData(normalized_eset), # This should be already correctly populated
      featureData = Biobase::AnnotatedDataFrame(final_feature_data_df),
      annotation = paste0(cdf_package_name, ".", cdf_version) # Store CDF info
    )

    # Verify sample names consistency (important)
    if(!identical(Biobase::sampleNames(final_eset), rownames(Biobase::pData(final_eset)))){
      warning("Sample name mismatch between expression data and phenoData for platform ", platform_id, ". This needs fixing.")
      # Attempt to fix if pheno_data_df has correct GSM IDs
      # This assumes order is preserved or GSM IDs in pheno_data_df$geo_accession match
      # sampleNames(final_eset)
      ordered_pd <- pheno_data_df[match(Biobase::sampleNames(final_eset), pheno_data_df$geo_accession),]
      rownames(ordered_pd) <- ordered_pd$geo_accession
      Biobase::phenoData(final_eset) <- Biobase::AnnotatedDataFrame(ordered_pd)
    }


    message("Successfully created ExpressionSet for platform ", platform_id)

    # --- Save RDS ---
    if (!is.null(output_dir) && !is.null(output_rds_filename)) {
      tryCatch({
        saveRDS(final_eset, file = output_rds_filename)
        message("Saved ExpressionSet to: ", output_rds_filename)
      }, error = function(e) {
        warning("Failed to save RDS file ", output_rds_filename, ": ", e$message)
      })
    }
    all_esets[[paste0(gse_accession, "_", platform_id)]] <- final_eset
  } # End loop over platforms in GSE

  message("\nProcessing of ", gse_accession, " complete.")
  if (length(all_esets) == 0) {
    message("No platforms were successfully processed for ", gse_accession)
    return(NULL)
  } else if (length(all_esets) == 1) {
    return(all_esets[[1]]) # Return single ExpressionSet if only one platform
  } else {
    return(all_esets) # Return a list of ExpressionSets if multiple platforms
  }
}
