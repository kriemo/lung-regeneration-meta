
## 2023-03-10
  - Reorganized oneDrive:
  .html reports are now in reports, with most recent results oin 2023-03-12 directory, prior versions in the "prior" directory
  the marker gene tables are placed in directories in results/2023-03-12/tables. 
  
  - plotted trajectories over clusters (See results/2023-03-12/05_figures.html)
  - recomputed markers using findMarkers from scran (See results/2023-03-12/tables/scanvi/clusters_9.xlsx)
  - plot umaps with max per study (See results/2023-03-12/05_figures.html)
  - started making publication figures (See results/2023-03-12/05_figures.html)
  - checked on cell #s from each study fall into cluster 7 from scvi and scanvi clusters (See results/2023-03-12/04_figures.html)
  
  - Scvi  Kobayashi_et_al_ncb (77) Riemondy_et_al_JCI-insight (112) Strunz_et_al_nc (59)
                  
  - scanvi Kobayashi_et_al_ncb (103) Riemondy_et_al_JCI-insight (219) Strunz_et_al_nc (244)
  
## 2023-02-24
  - Generated marker tables and figures comparing the transitional and aberrant basaloid populations
  - Built a UCSC cellbrowser for the scvi and scanvi integrations, hosted on aws S3 @  http://lung-meta-2023-zemans.s3-website-us-west-2.amazonaws.com/. 
  - Ran monocle3 as an additional pseuodotime algorithm to show.
  
## 2023-02-03
  - Performed pseudotime with slingshot on harmony, scanVI, and scVI projections. 

## 2023-01-14
  - Fixed heatmaps, was using an outdated singlecellexperiment object. 
  
## 2023-01-13
  - Made heatmaps of key genes from figure 7 for scANVI and scVI projections. 
  
## 2023-01-08
  - changed color palette to loupe-like colors
  - changed the coloring of the UMAP plots to plot in order of expression
  - looked into scVI AT1 cells, which are a mix of different clusters in each study,
    showing a discrepancy between the UMAP and clustering results.
