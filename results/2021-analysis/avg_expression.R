library(Seurat)
library(openxlsx)
so <- readRDS("data/Strunz_et_al_nc/so_whole.rds")
Idents(so) <- "cell.type"
so <- so[, !is.na(so$cell.type)]
struntz_whole_lung <- log1p(AverageExpression(so, return.seurat = FALSE)$RNA)

so <- readRDS("data/Strunz_et_al_nc/so_high_res.rds")
Idents(so) <- "cell_type_recombined"
so <- so[, !is.na(so$cell_type_recombined)]
struntz_high_res <- log1p(AverageExpression(so, return.seurat = FALSE)$RNA)

so <- readRDS("data/Habermann_et_al_bioRxiv/GSE135893_ILD_annotated_diet.rds")
DefaultAssay(so) <- "RNA"
so <- NormalizeData(so)
Idents(so) <- "celltype"
habermann <- log1p(AverageExpression(so, return.seurat = FALSE)$RNA)

res <- list(struntz_whole_lung = struntz_whole_lung,
            struntz_high_res = struntz_high_res,
            habermann = habermann)
res <- lapply(res, function(x) {
  x <- as.data.frame(x)
  x <- cbind(gene = rownames(x), x, row.names = NULL)
  class(x$gene) <- "Text"
  x
})
write.xlsx(res, "results/processed_data/average_expression_struntz_habermann_2022-08-19.xlsx")
