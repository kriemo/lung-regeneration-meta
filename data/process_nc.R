library(here)
library(data.table)
library(readxl)
library(Matrix)
library(stringr)
data_dir <- here("data")

#Struz et al NC.2020
# http://dx.doi.org/10.1038/s41467-020-17358-3

study_id <- "Strunz_et_al_nc"
dir.create(study_id, showWarnings = FALSE)
sce_fn <- file.path(data_dir, study_id, "sce.rds")

if(!file.exists(file.path(data_dir, study_id, "GSE141259_HighResolution_rawcounts.mtx.gz"))){
  options(timeout = 1e5)
  download.file("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE141nnn/GSE141259/suppl/GSE141259_HighResolution_rawcounts.mtx.gz",
                file.path(data_dir, study_id, "GSE141259_HighResolution_rawcounts.mtx.gz"))
  download.file("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE141nnn/GSE141259/suppl/GSE141259_HighResolution_genes.txt.gz",
                file.path(data_dir, study_id, "GSE141259_HighResolution_genes.txt.gz"))
  download.file("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE141nnn/GSE141259/suppl/GSE141259_HighResolution_barcodes.txt.gz",
                file.path(data_dir, study_id, "GSE141259_HighResolution_barcodes.txt.gz"))
  download.file("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE141nnn/GSE141259/suppl/GSE141259_HighResolution_cellinfo.csv.gz",
                file.path(data_dir, study_id, "GSE141259_HighResolution_cellinfo.csv.gz"))
}
mtx <- readMM(file.path(data_dir, study_id, "GSE141259_HighResolution_rawcounts.mtx.gz"))
colnames(mtx) <- readLines(file.path(data_dir, study_id, "GSE141259_HighResolution_genes.txt.gz"))
rownames(mtx) <- readLines(file.path(data_dir, study_id, "GSE141259_HighResolution_barcodes.txt.gz"))
mtx <- t(mtx)

mdata <- fread(file.path(data_dir, study_id, "GSE141259_HighResolution_cellinfo.csv.gz"),
               data.table = FALSE)
rownames(mdata) <- mdata$cell_barcode
stopifnot(all(rownames(mdata) %in% colnames(mtx)))
mtx <- mtx[, rownames(mdata)]
sce <- SingleCellExperiment(list(counts = mtx))
colData(sce) <- as(mdata, "DataFrame")
sce <- logNormCounts(sce)
saveRDS(sce, sce_fn)

