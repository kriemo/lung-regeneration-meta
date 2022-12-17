library(sceasy)
library(reticulate)
library(here)
library(SingleCellExperiment)
library(scater)
# generated a conda enviroment as follows
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

hvgs <- metadata(uc_sce)$hvgs$hvgs[1:2000]
sce <- uc_sce[hvgs, ]

adata <- sceasy::convertFormat(sce,
                               from="sce",
                               to="anndata",
                               outFile='hvg_uc_sce.h5ad',
                               main_layer = 'logcounts',
                               transfer_layers = 'counts'
                               )


# scVI --------------------------------

scvi$model$SCVI$setup_anndata(adata,
                              layer = 'counts',
                              batch_key = 'study')

# based on HCA lung atlas integration
model <- scvi$model$SCVI(adata, n_layers = 2L, n_latent = 20L, gene_likelihood="nb")
model$train()
latent = model$get_latent_representation()

scvi_mat <- as.matrix(latent)
rownames(scvi_mat) <- colnames(sce)
colnames(scvi_mat) <- paste0("LD", 1:ncol(scvi_mat))
reducedDim(uc_sce, "scvi") <- scvi_mat
set.seed(20221220)
scvi_sce <- runUMAP(uc_sce,
                    dimred = "scvi",
                    n_dimred = ncol(scvi_mat))

saveRDS(scvi_sce, here("results",
                       "2022-analysis",
                       "objects",
                       "scvi_sce.rds"))

# scANVI --------------------------------
lvae = scvi$model$SCANVI$from_scvi_model(
  model,
  adata=adata,
  unlabeled_category="Unknown",
  labels_key="coarse_cell_type"
)

# no label subsampling for this smaller dataset
lvae$train(max_epochs=20L)

latent_scanvi = lvae$get_latent_representation(adata)

scanvi_mat <- as.matrix(latent_scanvi)
rownames(scanvi_mat) <- colnames(uc_sce)
colnames(scanvi_mat) <- paste0("LD_AN", 1:ncol(scanvi_mat))

reducedDim(uc_sce, "scanvi") <- scanvi_mat

set.seed(20221220)
scanvi_sce <- runUMAP(uc_sce,
                    dimred = "scanvi",
                    n_dimred = ncol(scanvi_mat))

saveRDS(scanvi_sce, here("results",
                       "2022-analysis",
                       "objects",
                       "scanvi_sce.rds"))
