library(here)
library(data.table)
library(readxl)
library(Matrix)
library(stringr)
data_dir <- here("data")

# wu_et_al_cell
# http://dx.doi.org/10.1016/j.cell.2019.11.027

study_id <- "Wu_et_al_cell"
dir.create(study_id, showWarnings = FALSE)
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
mdata <- mdata[rownames(mdata) %in% colnames(mat), ]

mat <- mat[, rownames(mdata)]


  # also process at1 cells from previous pub (no metadata :( )
  at1_mat <- read_csv(file.path(data_dir,  study_id, "at1_data", "GSM2858341_AT1_P60.exprs.csv"))
  at1_mat <- column_to_rownames(at1_mat, "X1") %>%
    as.matrix(.) %>%
    as(., "sparseMatrix")
  so_at1 <- CreateSeuratObject(at1_mat)
  so_at1 <- FindVariableFeatures(so_at1) %>%
    ScaleData() %>%
    RunPCA(seed.value = 42) %>%
    RunUMAP(dims = 1:10, seed.value = 42) %>%
    FindNeighbors() %>%
    FindClusters(resolution = c(0.05, 0.1, 0.2), random.seed = 42)
  so_at1$cell_type <- ifelse(so_at1$RNA_snn_res.0.2 %in% c("0", "1", "2", "4"),
                             "AT1",
                             "other")
  so_at1 <- subset(so_at1,
                   subset = cell_type == "AT1")
  # get rid of a few outliers
  to_drop <- c("CAGGTGCCATTACGAC",
               "CTAGCCTGTGTTTGGT",
               "GTACGTAGTTCCGGCA",
               "TGGGAAGAGAAGAAGC")

  so_at1 <- subset(so_at1, cells = to_drop, invert = TRUE)

  so <- merge(so, so_at1)

  saveRDS(so, so_fn)
}
