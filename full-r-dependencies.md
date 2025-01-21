# R Dependencies

This project requires R version 4.3.2 or higher.

## Installation

You can install the required dependencies using the following R code:

```R
# Install if you don't have them
if (!require("BiocManager")) install.packages("BiocManager")
if (!require("remotes")) install.packages("remotes")

# Set Bioconductor version
BiocManager::install(version = "3.18")

### Install Bioconductor packages
bioc_packages <- c(
    "AnnotationDbi@1.64.1",
    "Biobase@2.62.0",
    "BiocGenerics@0.48.1",
    "BiocParallel@1.36.0",
    "BiocVersion@3.18.1",
    "DelayedArray@0.28.0",
    "DOSE@3.28.2",
    "EnhancedVolcano@1.20.0",
    "enrichplot@1.22.0",
    "fgsea@1.28.0",
    "GenomeInfoDbData@1.2.11",
    "GenomicRanges@1.54.1",
    "GO.db@3.18.0",
    "HDO.db@0.99.1",
    "IRanges@2.36.0",
    "KEGGREST@1.42.0",
    "limma@3.58.1",
    "MatrixGenerics@1.14.0",
    "qvalue@2.34.0",
    "S4Vectors@0.40.2",
    "SummarizedExperiment@1.32.0",
    "treeio@1.26.0",
    "XVector@0.42.0",
    "BiocFileCache@2.10.2",
    "biomaRt@2.58.2",
    "Biostrings@2.70.3",
    "clusterProfiler@4.10.1",
    "DESeq2@1.42.1",
    "edgeR@4.0.16",
    "GenomeInfoDb@1.38.8",
    "ggtree@3.10.1",
    "GOSemSim@2.28.1",
    "S4Arrays@1.2.1",
    "SparseArray@1.2.4",
    "zlibbioc@1.48.2"
)

### Install CRAN packages
cran_packages <- c(
    "abind@1.4-8",
    "ape@5.8-1",
    "aplot@0.2.4",
    "askpass@1.2.1",
    "backports@1.5.0",
    "bayestestR@0.15.0",
    "BH@1.87.0-1",
    "BiocManager@1.30.25",
    "bit@4.5.0.1",
    "bit64@4.5.2",
    "bitops@1.0-9",
    "blob@1.2.4",
    "boot@1.3-31",
    "broom@1.0.7",
    "cachem@1.1.0",
    "car@3.1-3",
    "carData@3.0-5",
    "cellranger@1.1.0",
    "cli@3.6.3",
    "codetools@0.2-20",
    "colorspace@2.1-1",
    "correlation@0.8.6",
    "corrplot@0.95",
    "cowplot@1.1.3",
    "cpp11@0.5.1",
    "crayon@1.5.3",
    "curl@6.1.0",
    "data.table@1.16.2",
    "datawizard@0.13.0",
    "DBI@1.2.3",
    "dbplyr@2.5.0",
    "Deriv@4.1.6",
    "digest@0.6.37",
    "doBy@4.6.24",
    "downloader@0.4",
    "dplyr@1.1.4",
    "effectsize@1.0.0",
    "EnvStats@3.0.0",
    "fansi@1.0.6",
    "farver@2.1.2",
    "fastmap@1.2.0",
    "fastmatch@1.1-6",
    "filelock@1.0.3",
    "formatR@1.14",
    "Formula@1.2-5",
    "fs@1.6.5",
    "futile.logger@1.4.3",
    "futile.options@1.0.1",
    "gargle@1.5.2",
    "generics@0.1.3",
    "ggforce@0.4.2",
    "ggfun@0.1.8",
    "ggh4x@0.3.0",
    "ggnewscale@0.5.0",
    "ggplot2@3.5.1",
    "ggplotify@0.1.2",
    "ggpubr@0.6.0",
    "ggraph@2.2.1",
    "ggrepel@0.9.6",
    "ggsci@3.2.0",
    "ggsignif@0.6.4",
    "glue@1.8.0",
    "googledrive@2.1.1",
    "googlesheets4@1.1.1",
    "graphlayouts@1.2.1",
    "gridExtra@2.3",
    "gridGraphics@0.5-1",
    "gson@0.1.0",
    "gtable@0.3.6",
    "hms@1.1.3",
    "httr@1.4.7",
    "ids@1.0.1",
    "igraph@2.1.3",
    "insight@1.0.0",
    "isoband@0.2.7",
    "jsonlite@1.8.9",
    "labeling@0.4.3",
    "lambda.r@1.2.4",
    "lattice@0.22-6",
    "lazyeval@0.2.2",
    "lifecycle@1.0.4",
    "lme4@1.1-35.5",
    "locfit@1.5-9.10",
    "magrittr@2.0.3",
    "MASS@7.3-60.0.1",
    "Matrix@1.6-5",
    "MatrixModels@0.5-3",
    "matrixStats@1.5.0",
    "memoise@2.0.1",
    "mgcv@1.9-1",
    "microbenchmark@1.5.0",
    "mime@0.12",
    "minqa@1.2.8",
    "modelbased@0.8.9",
    "modelr@0.1.11",
    "munsell@0.5.1",
    "nlme@3.1-166",
    "nloptr@2.1.1",
    "nnet@7.3-19",
    "nortest@1.0-4",
    "numDeriv@2016.8-1.1",
    "openssl@2.3.0",
    "parameters@0.24.0",
    "patchwork@1.3.0",
    "pbkrtest@0.5.3",
    "performance@0.12.4",
    "pillar@1.10.1",
    "pkgconfig@2.0.3",
    "plogr@0.2.0",
    "plyr@1.8.9",
    "png@0.1-8",
    "polyclip@1.10-7",
    "polynom@1.4-1",
    "pracma@2.4.4",
    "prettyunits@1.2.0",
    "progress@1.2.3",
    "purrr@1.0.2",
    "quantreg@5.99.1",
    "R6@2.5.1",
    "rappdirs@0.3.3",
    "RColorBrewer@1.1-3",
    "Rcpp@1.0.13-1",
    "RcppArmadillo@14.2.0-1",
    "RcppEigen@0.3.4.0.2",
    "RCurl@1.98-1.16",
    "readxl@1.4.3",
    "rematch@2.0.0",
    "rematch2@2.1.2",
    "renv@1.0.11",
    "reshape2@1.4.4",
    "rlang@1.1.4",
    "RSQLite@2.3.9",
    "rstatix@0.7.2",
    "scales@1.3.0",
    "scatterpie@0.2.4",
    "see@0.9.0",
    "shadowtext@0.1.4",
    "snow@0.4-4",
    "SparseM@1.84-2",
    "statmod@1.5.0",
    "stringi@1.8.4",
    "stringr@1.5.1",
    "survival@3.7-0",
    "svglite@2.1.3",
    "sys@3.4.3",
    "systemfonts@1.1.0",
    "tibble@3.2.1",
    "tidygraph@1.3.1",
    "tidyr@1.3.1",
    "tidyselect@1.2.1",
    "tidytree@0.4.6",
    "tweenr@2.0.3",
    "utf8@1.2.4",
    "uuid@1.2-1",
    "vctrs@0.6.5",
    "viridis@0.6.5",
    "viridisLite@0.4.2",
    "withr@3.0.2",
    "XML@3.99-0.18",
    "xml2@1.3.6",
    "yulab.utils@0.1.9"
)

# Install packages
BiocManager::install(bioc_packages)
install.packages(cran_packages)

# Install GitHub packages
remotes::install_github("wilkelab/ggridges@a8a99820")
```

## Alternative Installation Using renv

Alternatively, you can use `renv` for dependency management:

```R
# Install renv if not already installed
install.packages("renv")

# Initialize renv and restore dependencies
renv::init()
renv::restore()
```

## System Requirements

- R version 4.3.2 or higher
- Bioconductor 3.18
- Working internet connection for package downloads
- Appropriate system libraries for package compilation (varies by OS)

### Additional System Dependencies

Some packages may require system-level libraries. Here are common requirements by operating system:

#### Ubuntu/Debian
```bash
sudo apt-get update
sudo apt-get install build-essential libcurl4-openssl-dev libxml2-dev libssl-dev libfontconfig1-dev libharfbuzz-dev libfribidi-dev libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev
```

#### macOS (using Homebrew)
```bash
brew install openssl xml2 freetype libpng
```

#### Windows
- Rtools43 must be installed for your R version (4.3.x)
- Download from: https://cran.r-project.org/bin/windows/Rtools/

## Package Categories

### Core Analysis Packages
- DESeq2 (1.42.1) - Differential expression analysis
- edgeR (4.0.16) - Differential expression analysis
- limma (3.58.1) - Linear models for microarray/RNA-seq data
- clusterProfiler (4.10.1) - Enrichment analysis

### Visualization Packages
- ggplot2 (3.5.1) - Data visualization
- ggtree (3.10.1) - Phylogenetic tree visualization
- EnhancedVolcano (1.20.0) - Volcano plots
- enrichplot (1.22.0) - Enrichment visualization

### Data Manipulation
- dplyr (1.1.4) - Data manipulation
- tidyr (1.3.1) - Data tidying
- data.table (1.16.2) - Fast data manipulation

### Bioinformatics
- Biostrings (2.70.3) - Biological sequences manipulation
- GenomicRanges (1.54.1) - Genomic intervals manipulation
- biomaRt (2.58.2) - Access to biomart databases

## Notes
- All package versions are pinned to ensure reproducibility
- Some packages might have additional system dependencies depending on your OS
- It's recommended to use `renv` for reproducible environments
- The GitHub package `ggridges` is installed from a specific commit

