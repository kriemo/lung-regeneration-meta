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
    plotUMAP(sce[, sce$random_order],
             point_size = 0.5,
             colour_by = x,
             other_fields = "study") +
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
    plotUMAP(sce[, sce$random_order],
             point_size = 0.5,
             colour_by = x,
             other_fields = "study") +
      facet_wrap(~study)
  })
  names(ccaplots) <- cca_to_plot

  act_to_plot <- c("Lcn2", "Ly6i", "Il33")
  actplots <- lapply(act_to_plot, function(x){
    plotUMAP(sce[, sce$random_order],
             point_size = 0.5,
             colour_by = x,
             other_fields = "study") +
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
