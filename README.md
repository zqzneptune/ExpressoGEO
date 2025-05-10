# ExpressoGEO

`ExpressoGEO` is an R package for downloading, normalizing, and annotating Affymetrix microarray data from GEO.

## Installation

1. **Install `remotes`** :
    ```R
    if (!requireNamespace("remotes", quietly = TRUE))
        install.packages("remotes")
    ```
2.  **Install `BiocManager`** :
    ```R
    if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
    ```
3. **Install Bioconductor packages** :
    ```R
    BiocManager::install(c('GEOquery', 'Biobase', 'affyio', 'affy', 'oligo', 'AnnotationHub', 'ensembldb'))
    ```
    
4.  **Install `ExpressoGEO`**
    ```R
    remotes::install_github("zqzneptune/ExpressoGEO")
    ```
    *(This will also install required Bioconductor dependencies like `GEOquery`, `oligo`, `AnnotationHub`, etc.)*

## Quick Usage

```R
library(ExpressoGEO)

# Process a GEO Series (GSE). Replace "GSE12345" with a valid Affymetrix GSE.
# Results (ExpressionSet) will be saved in a temporary directory by default.
eset_result <- process_gse(gse_accession = "GSE12345")

# To specify an output directory:
# eset_result <- process_gse(gse_accession = "GSE12345", output_dir = "./my_geo_data")

if (!is.null(eset_result)) {
  print(eset_result)
}
