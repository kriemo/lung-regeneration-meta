library(here)
library(rmarkdown)
rd <- here("results")
rmds <- c("01_preprocess.Rmd",
          "02_coembed.Rmd",
          "03_pseudotime-slingshot.Rmd",
          "03_pseudotime-monocle3.Rmd",
          "04_additional_analysis.Rmd",
          "05_figures_scvi.Rmd")
rmd_paths <- file.path(rd, rmds)
lapply(rmd_paths[6], render)
