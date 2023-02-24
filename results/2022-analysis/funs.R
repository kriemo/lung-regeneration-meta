library(bluster)
get_correction_statistics <- function(sce,
                                      clusters,
                                      batch = "study",
                                      og_cells = "cell_type"){
  studies <- unique(sce$study)
  aris <- vector("list", 3)
  for(i in seq_along(studies)){
    s <- as.character(studies[i])
    aris[[i]] <- pairwiseRand(clusters[sce[[batch]] == s],
                              sce[[og_cells]][sce[[batch]] == s],
                              mode="index")
    names(aris)[i] <- s
  }
  aris
}


plot_correction_umaps <- function(sce){
  cloupe_cols <- rev(RColorBrewer::brewer.pal(11, "RdGy")[c(1:5, 7)])

  pstudy_split <- plotUMAP(sce,
                 point_size = 0.5,
                 colour_by = "cell_type",
                 other_fields = "study") +
    facet_wrap(~study)

  pstudy <- plotUMAP(sce,
                 point_size = 0.5,
                 colour_by = "study",
                 order_by = "random_order")

  p2 <- plotUMAP(sce,
                 point_size = 0.5,
                 colour_by = "cell_type",
                 order_by = "random_order")

  p3 <- plotUMAP(sce,
                 point_size = 0.5,
                 colour_by = "old_coarse_cell_types",
                 other_fields = "study") +
    facet_wrap(~study)

  to_plot <- c(
    "Sftpc",
    "Pdpn",
    "Scgb1a1",
    "Mki67",

    "Cldn4",
    "Sfn",
    "Hopx",
    "Igfbp2",
    "Krt8",

    "Foxj1"
  )
  gplots <- lapply(to_plot, function(x){
    idx <- order(logcounts(sce[x, ]))
    plotUMAP(sce[, idx],
             point_size = 0.5,
             colour_by = x,
             other_fields = "study") +
      scale_color_gradientn(name = x, colors = cloupe_cols) +
      facet_wrap(~study)
  })
  names(gplots) <- to_plot

  cca_to_plot <- c(
    "Cdkn1a",
    "Cdkn2a",
    "Cdkn2b",
    "Trp53",

    "Pdlim7",
    "Mcam",
    "Krt17",
    "Sox4",
    "Itga2",
    "Palld"
  )

  ccaplots <- lapply(cca_to_plot, function(x){
    idx <- order(logcounts(sce[x, ]))
    plotUMAP(sce[, idx],
             point_size = 0.5,
             colour_by = x,
             other_fields = "study") +
      scale_color_gradientn(name = x, colors = cloupe_cols) +
      facet_wrap(~study)
  })
  names(ccaplots) <- cca_to_plot

  act_to_plot <- c("Lcn2", "Ly6i", "Il33")
  actplots <- lapply(act_to_plot, function(x){
    idx <- order(logcounts(sce[x, ]))
    plotUMAP(sce[, idx],
             point_size = 0.5,
             colour_by = x,
             other_fields = "study") +
      scale_color_gradientn(name = x, colors = cloupe_cols) +
      facet_wrap(~study)
  })
  names(actplots) <- act_to_plot

  list(by_study = pstudy,
       by_study_split = pstudy_split,
       all = p2,
       all_old = p3,
       genes = gplots,
       cca_genes = ccaplots,
       act_genes = actplots)
}

# https://stackoverflow.com/a/54136863/6276041
quiet <- function(x) { 
  sink(tempfile()) 
  on.exit(sink()) 
  invisible(force(x)) 
} 

# from scater

sccol_pals <- list(
  tableau20 = c(
    "#1F77B4",
             "#AEC7E8",
             "#FF7F0E",
             "#FFBB78",
             "#2CA02C",
             "#98DF8A",
             "#D62728",
             "#FF9896",
             "#9467BD",
             "#C5B0D5",
             "#8C564B",
             "#C49C94",
             "#E377C2",
             "#F7B6D2",
             "#7F7F7F",
             "#C7C7C7",
             "#BCBD22",
             "#DBDB8D",
             "#17BECF",
             "#9EDAE5"
  ),
  tableau10medium = c(
    "#729ECE",
             "#FF9E4A",
             "#67BF5C",
             "#ED665D",
             "#AD8BC9",
             "#A8786E",
             "#ED97CA",
             "#A2A2A2",
             "#CDCC5D",
             "#6DCCDA"
  ),
  colorblind10 = c(
    "#006BA4",
             "#FF800E",
             "#ABABAB",
             "#595959",
             "#5F9ED1",
             "#C85200",
             "#898989",
             "#A2C8EC",
             "#FFBC79",
             "#CFCFCF"
  ),
  colourblind10 = c(
    "#006BA4",
             "#FF800E",
             "#ABABAB",
             "#595959",
             "#5F9ED1",
             "#C85200",
             "#898989",
             "#A2C8EC",
             "#FFBC79",
             "#CFCFCF"
  ),
  trafficlight = c(
    "#B10318",
             "#DBA13A",
             "#309343",
             "#D82526",
             "#FFC156",
             "#69B764",
             "#F26C64",
             "#FFDD71",
             "#9FCD99"
  ),
  purplegray12 = c(
    "#7B66D2",
             "#A699E8",
             "#DC5FBD",
             "#FFC0DA",
             "#5F5A41",
             "#B4B19B",
             "#995688",
             "#D898BA",
             "#AB6AD5",
             "#D098EE",
             "#8B7C6E",
             "#DBD4C5"
  ),
  bluered12 = c(
    "#2C69B0",
             "#B5C8E2",
             "#F02720",
             "#FFB6B0",
             "#AC613C",
             "#E9C39B",
             "#6BA3D6",
             "#B5DFFD",
             "#AC8763",
             "#DDC9B4",
             "#BD0A36",
             "#F4737A"
  ),
  greenorange12 = c(
    "#32A251",
             "#ACD98D",
             "#FF7F0F",
             "#FFB977",
             "#3CB7CC",
             "#98D9E4",
             "#B85A0D",
             "#FFD94A",
             "#39737C",
             "#86B4A9",
             "#82853B",
             "#CCC94D"
  ),
  cyclic = c(
    "#1F83B4",
             "#1696AC",
             "#18A188",
             "#29A03C",
             "#54A338",
             "#82A93F",
             "#ADB828",
             "#D8BD35",
             "#FFBD4C",
             "#FFB022",
             "#FF9C0E",
             "#FF810E",
             "#E75727",
             "#D23E4E",
             "#C94D8C",
             "#C04AA7",
             "#B446B3",
             "#9658B1",
             "#8061B4",
             "#6F63BB"
  )
)
