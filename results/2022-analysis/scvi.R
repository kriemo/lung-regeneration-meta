library(sceasy)
library(reticulate)
library(here)
library(SingleCellExperiment)
library(scater)

# generated a conda environment as follows
#
# mamba create -n scvi-env python=3.9
# conda activate scvi-env
# mamba install -c anaconda pip
# mamba install scvi-tools -c conda-forge
# mamba install -c conda-forge scanpy python-igraph leidenalg

# then installed sceasy
# devtools::install_github("cellgeni/sceasy")

sce_fn <-  here("results", "2022-analysis", "objects", "uc_sce.rds")
uc_sce <- readRDS(sce_fn)

use_condaenv('scvi-env')

sc <- import("scanpy", convert = FALSE)
scvi <- import("scvi", convert = FALSE)

hvgs <- metadata(uc_sce)$hvgs$hvgs[1:3000]
sce <- uc_sce[hvgs, ]

out_h5 <- here("results", "2022-analysis", "objects", "hvg_uc_sce.h5ad")
adata <- sceasy::convertFormat(sce,
                               from="sce",
                               to="anndata",
                               outFile=out_h5,
                               main_layer = 'logcounts',
                               transfer_layers = 'counts'
                               )

# use google colab for scvi, instead of via local R or python, as it provides access to GPU