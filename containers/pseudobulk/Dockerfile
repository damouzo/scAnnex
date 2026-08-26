# Pseudobulk DESeq2 container for scAnnex (DESeq2, apeglm, optparse, jsonlite)
# Based on Bioconductor 3.20 (R 4.4.x). Published to GHCR by build-containers.yml.

FROM bioconductor/bioconductor_docker:RELEASE_3_20

LABEL org.opencontainers.image.title="scannex-pseudobulk"
LABEL org.opencontainers.image.description="Pseudobulk DE for scAnnex: DESeq2, apeglm, optparse, jsonlite, ggplot2"
LABEL org.opencontainers.image.source="https://github.com/damouzo/scannex"
LABEL org.opencontainers.image.licenses="MIT"

# Bioconductor DE stack
RUN R -e " \
    BiocManager::install(c( \
        'DESeq2', \
        'apeglm' \
    ), ask=FALSE, update=FALSE, version=BiocManager::version()) \
"

# CRAN packages
RUN R -e " \
    install.packages(c( \
        'optparse', 'jsonlite', 'data.table', 'ggplot2' \
    ), repos='https://cloud.r-project.org', quiet=TRUE) \
"

# Sanity check
RUN Rscript -e " \
    suppressPackageStartupMessages({ \
        library(DESeq2); library(optparse); library(jsonlite); \
        library(data.table); library(ggplot2) \
    }); \
    cat('All pseudobulk packages loaded OK | DESeq2', as.character(packageVersion('DESeq2')), '\n') \
"