
## Notes on marker gene detection

To detect marker genes I've used an approach called findMarkers() from the scran bioconductor R package, which differs
from more commonly used approaches (e.g. those in the Seurat R pacakge). The common approach used to detect marker genes tests
for gene expression differences between a cluster vs all other cells in the dataset, this is then repeated for each cluster.
This approach can detect specific markers for a cluster, however it is sensitive to the composition of the sample. For example
if one cluster dominates the dataset, then most of the markers will be biased towards distinguishing between that cluster,
rather than between other clusters.

The findMarkers() approach conducts pair-wise comparisons between every cluster in the dataset, then provides a combined p-value that can be used to rank markers. This approach is more complicated to interpret, however allows tailoring the analysis
to the particular goals in mind. For example, we can identify genes that globally distinguish one cluster from all other clusters
(very specific marker genes), or we can identify genes that distinguish at least 1 cluster from another cluster (very sensitive
to markers between subpopulations of cells), or a combination of the two approaches, e.g. requiring that a marker gene is differentially expressed in at least 50% of pair-wise comparisions.

For this analysis I have ranked genes based on `FDR`, followed by `summary.logFC` in the case of ties. The analysis has been done in a manner that requires 1/3 of pairwise comparisons to be significant in order to be reported (`pval.type = "some"`).

# Column definitions for marker gene files

`summary.logFC`: log2FC of the third most significant cluster comparisons
`mean.logFC`: mean log2FC of the pairwise cluster comparisons

`self.average`: The average log-normalized abundance of a gene in the indicated cluster
`other.average`: The average log-normalized abundance of a gene in all other cells
`self.detected`: The proportion of cells in the indicated cluster expressing the gene
`other.detected`: The proportion of cells not in the cluster expressing the gene.


Material quoted or reinterpreted from the following references:

Amezquita, R.A., Lun, A.T.L., Becht, E. et al. Orchestrating single-cell analysis with Bioconductor. Nat Methods 17, 137–145 (2020). https://doi.org/10.1038/s41592-019-0654-x

Lun ATL, McCarthy DJ, Marioni JC (2016). “A step-by-step workflow for low-level analysis of single-cell RNA-seq data with Bioconductor.” _F1000Res._, *5*, 2122. doi:10.12688/f1000research.9501.2 <https://doi.org/10.12688/f1000research.9501.2>.

The scran Bioconductor package documentation:
https://bioconductor.org/packages/devel/bioc/html/scran.html