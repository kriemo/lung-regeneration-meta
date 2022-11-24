library(here)
library(data.table)
library(readxl)
library(Matrix)
library(stringr)
data_dir <- here("data")

# our jci data
# https://doi.org/10.1172/jci.insight.123637

study_id <- "Riemondy_et_al_JCI-insight"
dir.create(study_id, showWarnings = FALSE)
sce_fn <- file.path(data_dir, study_id, "sce.rds")

if(!file.exists(file.path(data_dir, study_id, "GSE113049_count_matrix.tsv.gz"))){
  download.file("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE113nnn/GSE113049/suppl/GSE113049_count_matrix.tsv.gz",
                file.path(data_dir, study_id, "GSE113049_count_matrix.tsv.gz"))
  download.file("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE113nnn/GSE113049/suppl/GSE113049_cell_metadata.tsv.gz",
                file.path(data_dir, study_id, "GSE113049_cell_metadata.tsv.gz"))
}
mat <- fread(file.path(data_dir, study_id, "GSE113049_count_matrix.tsv.gz"),
             data.table = FALSE)
rownames(mat) <- mat[, 1]
mat[, 1] <- NULL
mat <- as.matrix(mat) |> as("sparseMatrix")

mdata <- fread(file.path(data_dir, study_id, "GSE113049_cell_metadata.tsv.gz"),
               data.table = FALSE)
rownames(mdata) <- mdata$cell
cell_id_relabel <- c("Endothelial/Fibroblast" = "Endothelial/Fibroblast",
                     "Basal" = "Basal",
                     "Club" = "Club",
                     "Ciliated" = "Ciliated",
                     "Naive AEC1" = "Naive AEC1",
                     "Macrophage" = "Macrophage",
                     "Other Injured AEC2" = "Other Injured AEC2",
                     "Injured AEC2: Proliferating" = "Injured AEC2: Proliferating",
                     "Injured AEC2: Cell Cycle Arrest" = "Injured AEC2: Early Transdifferentiating",
                     "Injured AEC2: Transdifferentiating" = "Injured AEC2: Late Transdifferentiating",
                     "Naive AEC2" = "Naive AEC2")

stopifnot(all(mdata$cell_type %in% names(cell_id_relabel)))

mdata$pretty_cell_labels <- factor(cell_id_relabel[mdata$cell_type], levels = cell_id_relabel)

stopifnot(all(rownames(mdata) %in% colnames(mat)))
mat <- mat[, rownames(mdata)]
sce <- SingleCellExperiment(list(counts = mat))
colData(sce) <- as(mdata, "DataFrame")
sce <- logNormCounts(sce)
saveRDS(sce, sce_fn)

