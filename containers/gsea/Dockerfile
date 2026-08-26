# GSEA container for scAnnex (clusterProfiler, ReactomePA, org.Hs.eg.db, org.Mm.eg.db)
# Based on Bioconductor 3.20 (R 4.4.x). Published to GHCR by build-containers.yml.

FROM bioconductor/bioconductor_docker:RELEASE_3_20

LABEL org.opencontainers.image.title="scannex-gsea"
LABEL org.opencontainers.image.description="GSEA analysis for scAnnex: clusterProfiler, ReactomePA, enrichplot, org.Hs.eg.db, org.Mm.eg.db, pathview"
LABEL org.opencontainers.image.source="https://github.com/damouzo/scannex"
LABEL org.opencontainers.image.licenses="MIT"

# Bioconductor annotation DBs + enrichment stack (full, non-pruned versions)
RUN R -e " \
    BiocManager::install(c( \
        'clusterProfiler', \
        'enrichplot', \
        'DOSE', \
        'ReactomePA', \
        'pathview', \
        'KEGGREST', \
        'org.Hs.eg.db', \
        'org.Mm.eg.db', \
        'GO.db', \
        'AnnotationDbi', \
        'fgsea' \
    ), ask=FALSE, update=FALSE, version=BiocManager::version()) \
"

# CRAN packages for plotting + IO
RUN R -e " \
    install.packages(c( \
        'data.table', 'ggplot2', 'dplyr', 'stringr', 'jsonlite', 'ggridges' \
    ), repos='https://cloud.r-project.org', quiet=TRUE) \
"

# Sanity check
RUN Rscript -e " \
    suppressPackageStartupMessages({ \
        library(clusterProfiler); library(enrichplot); library(ggridges); \
        library(ReactomePA); library(org.Hs.eg.db); library(org.Mm.eg.db); \
        library(pathview); library(ggplot2); library(data.table) \
    }); \
    cat('All GSEA packages loaded OK\n') \
"