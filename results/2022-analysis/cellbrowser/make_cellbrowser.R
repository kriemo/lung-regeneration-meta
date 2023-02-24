library(Seurat)
library(SingleCellExperiment)
library(scbp)
library(here)
library(data.table)
library(readr)
scvi_sce <- readRDS(here("results", "2022-analysis", "objects", "20221230_scvi_sce.rds"))
scanvi_sce <- readRDS(here("results", "2022-analysis", "objects", "20221230_scanvi_sce.rds"))


write_config <- function(x, path){
  paste0(names(x), '="', x, '"') |>
    writeLines(path)
}

cb_outdir <- here("results", "2022-analysis", "cellbrowser")

collection_cb <- c(shortLabel="Meta-analysis of murine lung regeneration")
collection_desc <- c(title="Meta-analysis of murine lung regeneration",
                     abstract="10x Genomics libraries were processed using Seurat and Bioconductor packages, and sc-vi tools to generate unified UMAPs and clusters",
                     unitDesc="Log-normalized counts")
write_config(collection_cb, file.path(cb_outdir, "cellbrowser.conf"))
write_config(collection_desc, file.path(cb_outdir, "desc.conf"))
# Sys.setenv("CBOUT" = file.path(cb_outdir,"lung-meta-cb"))

writeLines(paste0("dataRoot=", sQuote(cb_outdir)), "~/.cellbrowser.conf")

prep_cb <- function(sce, 
                    cluster_col = "clusters_8", 
                    study_col = "study", 
                    integration = "scvi",
                    mrk_file = here("results/2022-analysis/tables/scvi/clusters_8_markers.csv")){
  so <- as.Seurat(sce)
  so[["RNA"]] <- so[["originalexp"]]
  cids <- split(colnames(so), so$study)
  embeddings <- c("UMAP")
  for(i in seq_along(cids)) {
    study <- names(cids)[i]
    study_cells <- cids[[i]]
    rid <- paste0("UMAP_", study)
    embeddings <- append(embeddings, rid)
    so@reductions[[rid]] <- so@reductions$UMAP
    Key(so@reductions[[rid]]) <- paste0(rid, "_")
    
    to_drop <- colnames(so)[!colnames(so) %in% study_cells]
    so@reductions[[rid]]@cell.embeddings[to_drop, ] <- NA
    
  }

  cols_to_keep <- c(
    `genes per cell` = "nFeature_originalexp",
    `UMIs per cell` = "nCount_originalexp",
    `Published cell type` = "cell_type",
    `Study` = "study",
    `Clusters` = cluster_col)
  
  if(endsWith(mrk_file, ".csv")){
    tf <- tempfile(fileext = ".tsv")
    head(fread(mrk_file), 100) |>
      fwrite(tf, 
             quote = FALSE,
             sep = "\t",
             col.names = TRUE,
             row.names = FALSE)
    on.exit(unlink(tf))
    mrk_file <- tf
  }
  make_cellbrowser(so,
                   column_list = cols_to_keep,
                   secondary_cols = c("Study"),
                   project = paste0("lung-meta-", integration),
                   outdir = cb_outdir,
                   marker_file = mrk_file,
                   ident = "Clusters",
                   embeddings = embeddings,
                   skip_expr_matrix = TRUE,
                  # config = list(priority = 1),
                   description = list(
                     title = paste0(integration, " integration"),
                     description = paste0("10x genomics data integrated with ", integration)
                   ),
                   cellbrowser_dir = "/usr/local/bin/"
  )
}

# scvi
prep_cb(scvi_sce, 
        cluster_col = "clusters_8",
        study_col = "study", 
        integration = "scvi", 
        mrk_file = here("results/2022-analysis/tables/scvi/clusters_8_markers.csv"))

# scanvi
prep_cb(scanvi_sce, 
        cluster_col = "clusters_9",
        study_col = "study", 
        integration = "scanvi", 
        mrk_file = here("results/2022-analysis/tables/scanvi/clusters_9_markers.csv"))



build_cellbrowser(c(file.path(cb_outdir, "lung-meta-scvi/cellbrowser.conf"),
                    file.path(cb_outdir, "lung-meta-scanvi/cellbrowser.conf")),
                  file.path(cb_outdir,"lung-meta-cb"),
                  "/usr/local/bin/cbBuild")
