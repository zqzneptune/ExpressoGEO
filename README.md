# ExpressoGEO: Streamlined Affymetrix GEO Data Processing

`ExpressoGEO` is an R package designed to simplify and automate the process of downloading, normalizing, and annotating Affymetrix microarray data from the Gene Expression Omnibus (GEO). It handles the complexities of fetching raw CEL files, installing necessary custom Chip Definition File (CDF) packages, performing Robust Multi-array Average (RMA) normalization, and mapping probe IDs to Ensembl gene IDs and gene names.

## Features

*   **Automated GEO Download:** Fetches GEO Series (GSE) metadata and associated raw (`_RAW.tar`) CEL files.
*   **Custom CDF Management:** Automatically attempts to download and install appropriate Ensembl gene-based custom CDF packages from [mbni.org/customcdf](http://mbni.org/customcdf/) for supported Affymetrix platforms.
*   **RMA Normalization:** Performs RMA normalization using the `oligo` package.
*   **Gene Annotation:** Annotates probe sets to Ensembl gene IDs and gene names using `AnnotationHub` and `ensembldb`.
*   **ExpressionSet Output:** Returns processed data as a standard Bioconductor `ExpressionSet` object, ready for downstream analysis.
*   **Handles Multi-Platform GSEs:** Can process GSE series containing data from multiple microarray platforms, returning a list of `ExpressionSet` objects.

## Installation

To use `ExpressoGEO`, you need to have R (>= 4.0) installed. This package relies on several Bioconductor packages. The recommended way to install `ExpressoGEO` and its dependencies is via `BiocManager`.

1.  **Install `BiocManager`** (if you haven't already):
    ```R
    if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
    ```

2.  **Install `ExpressoGEO`**:

    *   **From GitHub (Recommended for the latest development version):**
        Replace `YourGitHubUsername` with the actual GitHub username/organization where the package is hosted.
        ```R
        BiocManager::install("YourGitHubUsername/ExpressoGEO")
        ```
    *   **From CRAN (Once submitted and accepted):**
        ```R
        # BiocManager::install("ExpressoGEO")
        ```
    *   **From a local source file (`.tar.gz`):**
        ```R
        # BiocManager::install("path/to/ExpressoGEO_version.tar.gz")
        ```

`BiocManager` will handle the installation of all required CRAN and Bioconductor dependencies, such as `GEOquery`, `Biobase`, `oligo`, `AnnotationHub`, `ensembldb`, etc.

## Usage Example

The primary function in this package is `process_gse()`.

```R
# Load the package
library(ExpressoGEO)

# Specify the GSE accession number you want to process
# IMPORTANT: Choose a GSE that uses Affymetrix platforms likely supported by
# mbni.org custom CDFs (e.g., MoGene, HuGene, Clariom S/D arrays).
# For this example, let's assume "GSE12345" is such a dataset.
# Replace "GSE12345" with a real GSE accession for testing.
gse_id <- "GSE12345" # <--- !!! REPLACE WITH A VALID GSE FOR YOUR TEST !!!

# Specify an output directory for downloaded files and the final RDS object
# Using a temporary directory for this example:
output_path <- tempfile(pattern = "ExpressoGEO_output_")
dir.create(output_path)
message("Output will be saved in: ", output_path)

# --- Process the GSE data ---
# This can take some time depending on the size of the GSE and internet speed.
# It will download data, potentially install CDFs, normalize, and annotate.
eset_list <- process_gse(
  gse_accession = gse_id,
  output_dir = output_path,
  force_reprocess = FALSE # Set to TRUE to re-process even if RDS exists
)

# --- Check the results ---
if (!is.null(eset_list)) {
  if (is(eset_list, "ExpressionSet")) {
    # Single platform processed
    message("Successfully processed one platform for ", gse_id)
    print(eset_list)
    # Explore the ExpressionSet
    # Biobase::exprs(eset_list)[1:5, 1:3] # View first few rows/cols of expression data
    # Biobase::pData(eset_list)[1:3, ]    # View first few rows of phenotype data
    # Biobase::fData(eset_list)[1:5, ]    # View first few rows of feature data
  } else if (is.list(eset_list) && all(sapply(eset_list, is, "ExpressionSet"))) {
    # Multiple platforms processed
    message("Successfully processed ", length(eset_list), " platforms for ", gse_id)
    for (platform_name in names(eset_list)) {
      message("\n--- Platform: ", platform_name, " ---")
      print(eset_list[[platform_name]])
    }
  }
  message("Processed data (RDS file) saved in: ", output_path)
} else {
  message("Processing failed or no suitable platforms found for ", gse_id)
  message("Check console messages for details. Common reasons include:")
  message("- Unsupported microarray platform (no custom CDF available).")
  message("- Issues downloading data from GEO or custom CDFs.")
  message("- Inconsistent metadata in the GEO record (e.g., mixed species/platforms).")
}

# Clean up the temporary directory (optional)
# unlink(output_path, recursive = TRUE, force = TRUE)
```

### Function Arguments

The `process_gse` function has several key arguments:

*   `gse_accession`: (Character) The GSE accession number (e.g., `"GSE65656"`). **Required.**
*   `output_dir`: (Character) Directory to save downloaded files and the final `ExpressionSet` (as an RDS file). Defaults to `tempdir()`. If `NULL`, the `ExpressionSet` is returned but not saved to disk.
*   `cdf_version`: (Character) Version of the custom CDF to use from mbni.org. Defaults to `"25.0.0"`.
*   `cdf_gene_map`: (Character) Gene mapping type for the custom CDF (e.g., `"ensg"` for Ensembl gene IDs). Defaults to `"ensg"`.
*   `force_reprocess`: (Logical) If `TRUE`, re-processes the data even if an RDS file already exists in `output_dir`. Defaults to `FALSE`.

## Expected Results

If processing is successful, the `process_gse()` function will:

1.  **Return Value:**
    *   An `ExpressionSet` object if only one microarray platform was processed from the GSE.
    *   A named `list` of `ExpressionSet` objects if the GSE contained multiple processable platforms. Each element in the list corresponds to a platform.
    *   `NULL` if processing failed for all platforms within the GSE.

2.  **Output File (if `output_dir` is specified):**
    *   An RDS file (e.g., `mm_MoGene-1_0-st-v1_GSE12345_GPL6246.RDS`) will be saved in the specified `output_dir` (or a subdirectory named after the GSE within `output_dir`). The filename typically includes species code, platform name, GSE ID, and GPL ID. This file contains the serialized `ExpressionSet` object.

### Contents of the `ExpressionSet`

Each `ExpressionSet` object will contain:

*   **Assay Data (`exprs(eset)`):** A matrix of normalized expression values (log2 scale after RMA). Rows typically correspond to Ensembl gene IDs (or probe IDs if annotation fails), and columns correspond to samples (GSM IDs).
*   **Pheno Data (`pData(eset)`):** A data frame containing sample metadata extracted from GEO (e.g., `geo_accession`, `title`, `source_name_ch1`, `organism_ch1`, characteristics).
*   **Feature Data (`fData(eset)`):** A data frame containing annotations for the features (rows of the expression matrix). Minimally, it will have `gene_id`. If annotation is successful, it will also include `gene_name` (e.g., HUGO symbols for human). Rownames will be the `gene_id`.
*   **Annotation (`annotation(eset)`):** A string indicating the custom CDF package name and version used for processing (e.g., `"pdrg230pmrnensg.25.0.0"`).

## Important Considerations & Troubleshooting

*   **Internet Connection:** A stable internet connection is required to download data from GEO, AnnotationHub, and custom CDF repositories.
*   **Platform Support:** The automated custom CDF installation relies on the availability of Ensembl gene-based CDFs at [mbni.org/customcdf/](http://mbni.org/customcdf/). Not all Affymetrix platforms are supported there. If a platform is not supported, processing for that platform will be skipped. Check the console output for messages regarding CDF availability.
*   **Disk Space:** Raw GEO data can be large. Ensure you have sufficient disk space in the `output_dir`.
*   **Processing Time:** Processing can be time-consuming, especially for large GSEs or when downloading/installing CDFs for the first time.
*   **Mixed GSEs:** If a GSE series contains samples from different species or vastly different microarray types on the same declared platform, `ExpressoGEO` will attempt to process homogeneous subsets but may skip heterogeneous ones.
*   **Firewall/Proxy Issues:** If you are behind a strict firewall or proxy, downloads from FTP (GEO) or HTTP (mbni.org, AnnotationHub) might fail. Configure your R environment for proxy use if necessary.
*   **AnnotationHub Cache:** `AnnotationHub` caches data locally. The first time you run it for a specific organism/database, it might take longer to download the annotation database.

## Related Information

*   **GEO (Gene Expression Omnibus):** [www.ncbi.nlm.nih.gov/geo/](https://www.ncbi.nlm.nih.gov/geo/)
*   **Bioconductor:** [www.bioconductor.org](https://www.bioconductor.org)
    *   `GEOquery` package
    *   `Biobase` package (for `ExpressionSet`)
    *   `oligo` package
    *   `AnnotationHub` and `ensembldb` packages
*   **Custom CDFs (BrainArray):** [Download Page](http://brainarray.mbni.med.umich.edu/Brainarray/Database/CustomCDF/CDF_download.asp)

## Contributing

[Optional: Add guidelines if you want others to contribute, e.g., bug reports, feature requests, pull requests.]
Issues and pull requests are welcome. Please check the GitHub repository's issues page.

## License

This package is licensed under the GPL-3 License. See the `LICENSE` file for details.
