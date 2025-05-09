# R/utils_annotation.R

#' Fetch Gene Annotations
#'
#' Fetches gene annotations (gene_id, gene_name) using AnnotationHub.
#'
#' @param species_code Two-letter species code (e.g., "hs", "mm", "rn").
#' @param gene_ids A character vector of gene IDs (e.g., Ensembl IDs).
#' @return A data frame with `gene_id` and `gene_name` columns.
#' @importFrom AnnotationHub AnnotationHub
#' @importFrom ensembldb genes
#' @importFrom dplyr left_join
#' @keywords internal
fetch_gene_annotations <- function(species_code, gene_ids) {
  # These could be moved to zzz.R or passed as arguments if more flexibility is needed
  ah_map <- c("mm" = "AH89211", "hs" = "AH89180", "rn" = "AH89282")

  if (!species_code %in% names(ah_map)) {
    stop("Unsupported species code for AnnotationHub: ", species_code)
  }

  message("Fetching annotations for species: ", species_code, " using AnnotationHub ID: ", ah_map[[species_code]])
  ah <- AnnotationHub::AnnotationHub()
  edb <- ah[[ah_map[[species_code]]]]

  # Check if gene_id and gene_name columns exist
  # This can vary slightly between EnsDb versions/objects
  all_genes_data <- tryCatch({
    ensembldb::genes(edb, columns = c("gene_id", "gene_name"), return.type = "data.frame")
  }, error = function(e) {
    # Fallback if 'gene_name' isn't directly available or other common names
    message("Primary columns 'gene_id', 'gene_name' not found or error. Trying 'symbol' for gene_name.")
    # You might need to inspect edb object to find appropriate columns
    # This is a common variation:
    alt_cols <- ensembldb::listColumns(edb, "gene")
    name_col <- if ("gene_name" %in% alt_cols) "gene_name" else if ("symbol" %in% alt_cols) "symbol" else NULL
    if(is.null(name_col)) stop("Could not find a suitable gene name column in EnsDb object.")

    ensembldb::genes(edb, columns = c("gene_id", name_col), return.type = "data.frame")
  })

  # Ensure column names are consistent for merging
  if ("symbol" %in% colnames(all_genes_data) && !"gene_name" %in% colnames(all_genes_data)) {
    colnames(all_genes_data)[colnames(all_genes_data) == "symbol"] <- "gene_name"
  }
  if (!all(c("gene_id", "gene_name") %in% colnames(all_genes_data))) {
    stop("EnsDb object does not contain expected 'gene_id' and 'gene_name' (or 'symbol') columns.")
  }

  # Create a data frame from input gene_ids for joining
  feature_df <- data.frame(gene_id = as.character(gene_ids), stringsAsFactors = FALSE)

  annotated_features <- dplyr::left_join(feature_df, all_genes_data[, c("gene_id", "gene_name")], by = "gene_id")

  # For genes not found in EnsDb, gene_name will be NA. Keep them.
  # Ensure rownames are gene_ids for ExpressionSet featureData
  rownames(annotated_features) <- annotated_features$gene_id
  return(annotated_features)
}
