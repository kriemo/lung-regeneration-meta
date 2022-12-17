library(here)
library(data.table)
library(readxl)
library(Matrix)
library(stringr)
library(SingleCellExperiment)
library(scran)
library(scater)
library(BiocSingular)

data_dir <- here("data")

# wu_et_al_cell
# http://dx.doi.org/10.1016/j.cell.2019.11.027

study_id <- "Wu_et_al_cell"
dir.create(file.path(data_dir, study_id), showWarnings = FALSE)
sce_fn <- file.path(data_dir, study_id, "sce.rds")

mat_fn <- file.path(data_dir, study_id, "GSE138585_C0_C21_N0_N21.UMI.csv.gz")
if(!file.exists(mat_fn)){
  download.file("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE138nnn/GSE138585/suppl/GSE138585_C0_C21_N0_N21.UMI.csv.gz",
                file.path(data_dir, study_id, "GSE138585_C0_C21_N0_N21.UMI.csv.gz"))
}
mat <- fread(file.path(data_dir,  study_id, "GSE138585_C0_C21_N0_N21.UMI.csv.gz"), data.table = FALSE)
rownames(mat) <- mat[, 1]
mat[, 1] <- NULL
mat <- as.matrix(mat) |> as("sparseMatrix")

# provided by the authors
mdata <- fread(file.path(data_dir, study_id, "meta_info_C0_C7_C21.csv"), data.table = FALSE)
rownames(mdata) <- mdata[, 1]
mdata[, 1] <- NULL
mdata <- mdata[rownames(mdata) %in% colnames(mat), c("sample_name", "cell_type")]
mat <- mat[, rownames(mdata)]

sce <- SingleCellExperiment(list(counts = mat))
sce <- logNormCounts(sce)
colData(sce) <- cbind(colData(sce), mdata)

# also process at1 cells from previous pub (no metadata)
mat_fn <- file.path(data_dir, study_id, "GSM2858341_AT1_P60.exprs.csv.gz")

if(!file.exists(mat_fn)){
  download.file( "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2858nnn/GSM2858341/suppl/GSM2858341_AT1_P60.exprs.csv.gz",
                mat_fn)
}

at1_mat <- fread(mat_fn, data.table = FALSE)
rownames(at1_mat) <- at1_mat[, 1]
at1_mat[, 1] <- NULL
at1_mat <- as.matrix(at1_mat) |> as("sparseMatrix")

seed_val <- 20221125
sce_at1 <- SingleCellExperiment(list(counts = at1_mat))
sce_at1 <- logNormCounts(sce_at1)
dec <- modelGeneVarByPoisson(sce_at1)
top_feats <- getTopHVGs(sce_at1, prop=0.1)

set.seed(seed_val)
sce_at1 <- denoisePCA(sce_at1, technical=dec, subset.row=top_feats)
set.seed(seed_val)
sce_at1 <- runUMAP(sce_at1, dimred="PCA")
set.seed(seed_val)
snn.gr <- buildSNNGraph(sce_at1, use.dimred="PCA", k=20)
set.seed(seed_val)
colLabels(sce_at1) <- factor(igraph::cluster_louvain(snn.gr, resolution = 0.1)$membership)

sce_at1$cell_type <- ifelse(sce_at1$label == "1",
                           "AT1",
                           "other")
sce_at1 <- sce_at1[, sce_at1$cell_type == "AT1"]

reducedDims(sce_at1) <- NULL
sce_at1$sample_name <- "d60"
sce_at1$label <- NULL

genes_shared <- intersect(rownames(sce), rownames(sce_at1))
sce <- cbind(sce[genes_shared, ], sce_at1[genes_shared, ])

saveRDS(sce, sce_fn)

