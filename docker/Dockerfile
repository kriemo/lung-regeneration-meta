FROM ghcr.io/bioconductor/bioconductor:release_3_16
ENV BIOCONDUCTOR_USE_CONTAINER_REPOSITORY=FALSE
WORKDIR /home/rstudio
COPY r-pkgs.R r-pkgs.R
RUN R -e 'source("r-pkgs.R")'