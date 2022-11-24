library(here)
library(data.table)
library(readxl)
library(Matrix)
library(stringr)
data_dir <- here("data")

# Kobayashi_et_al_ncb
# https://doi.org/10.1038/s41556-020-0542-8

study_id <- "Kobayashi_et_al_ncb"
dir.create(study_id, showWarnings = FALSE)
sce_fn <- file.path(data_dir, study_id, "sce.rds")
mat_fn <- file.path(data_dir, study_id, "GSM4210295_MTECplus.tsv.gz")
if(!file.exists(mat_fn)){
  download.file("https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4210nnn/GSM4210295/suppl/GSM4210295_MTECplus.tsv.gz",
                mat_fn)
}

mat <- fread(mat_fn, data.table = FALSE)
rownames(mat) <- mat[, 1]
mat[, 1] <- NULL
mat <- as.matrix(mat) |> as("sparseMatrix")

# excel file provided by authors
mdata <- read_excel(file.path(data_dir,
                              study_id,
                              "epithelial cells_MTEC_organoids_d10_active_ident.xlsx"))
mdata$cell_id <- str_split(mdata$Cell_barcode,
                           ":",
                           simplify = TRUE)[, 2] |> str_remove("x$")
mdata <- as.data.frame(mdata)
rownames(mdata) <- mdata$cell_id
mdata <- mdata[, "active.ident", drop = FALSE]
colnames(mdata) <- "cell_type"

stopifnot(all(rownames(mdata) %in% colnames(mat)))
mat <- mat[, rownames(mdata)]
sce <- SingleCellExperiment(list(counts = mat))
colData(sce)$cell_type <- mdata$cell_type
sce <- logNormCounts(sce)
saveRDS(sce, sce_fn)
