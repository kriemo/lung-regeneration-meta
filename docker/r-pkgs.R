BiocManager::install(version = "3.16", ask = FALSE)
pkgs <- c("openxlsx",
          "Seurat",
          "rmarkdown",
          "knitr",
          "ComplexHeatmap",
          "cowplot",
          "dplyr",
          "here",
          "readxl",
          "tidyverse",
          "biomaRt",
          "tibble",
          "Matrix",
          "eulerr",
          "viridis",
          "ggrepel",
          "data.table",
          "scater",
          "scran",
          "scuttle",
          "SingleCellExperiment",
          "batchelor",
          "BiocParallel",
          "stringr",
          "igraph",
          "purrr",
          "slingshot",
          "clustree",
          "bluster",
          "RColorBrewer",
          "magrittr",
          "reticulate",
          "clustifyr",
          "immunogenomics/harmony",
          "immunogenomics/presto",
          "cellgeni/sceasy",
          "rnabioco/scbp")

BiocManager::install(pkgs, Ncpus = 6, version = "3.16")

monocle3_pkgs <- c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
                                       'limma', 'lme4', 'S4Vectors', 'SingleCellExperiment',
                                       'SummarizedExperiment', 'batchelor', 'HDF5Array',
                                       'terra', 'ggrastr', 'grr')
BiocManager::install(monocle3_pkgs, Ncpus = 6, version = "3.16")

# need to also use archived version of Matrix.utils, as it is no longer on CRAN
install.packages("https://cran.r-project.org/src/contrib/Archive/Matrix.utils/Matrix.utils_0.9.8.tar.gz", dependencies = TRUE, type = "source", repos = NULL)

BiocManager::install(c('cole-trapnell-lab/monocle3', 
                       'satijalab/seurat-wrappers'),
                     Ncpus = 6, version = "3.16")


