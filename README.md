# ExpressGEO

`ExpressGEO` is an R package for downloading, normalizing, and annotating Affymetrix microarray data from GEO.

## Installation

1.  **Install `BiocManager`** (if you don't have it):
    ```R
    if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
    ```

2.  **Install `ExpressGEO`**:
    (Replace `YourGitHubUsername` with the actual GitHub username/organization)
    ```R
    BiocManager::install("zqzneptune/ExpressGEO")
    ```
    *(This will also install required Bioconductor dependencies like `GEOquery`, `oligo`, `AnnotationHub`, etc.)*

## Quick Usage

```R
library(ExpressGEO)

# Process a GEO Series (GSE). Replace "GSE12345" with a valid Affymetrix GSE.
# Results (ExpressionSet) will be saved in a temporary directory by default.
eset_result <- process_gse(gse_accession = "GSE12345")

# To specify an output directory:
# eset_result <- process_gse(gse_accession = "GSE12345", output_dir = "./my_geo_data")

if (!is.null(eset_result)) {
  print(eset_result)
}
