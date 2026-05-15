FROM bioconductor/bioconductor_docker:RELEASE_3_20

LABEL org.opencontainers.image.title="scannex-dashboard"
LABEL org.opencontainers.image.description="R/Python environment for scAnnex scRNA-seq interactive dashboard (Shiny + scanpy/AnnData)"
LABEL org.opencontainers.image.version="1.0.0"
LABEL org.opencontainers.image.licenses="MIT"

# ---- System dependencies -------------------------------------------------------
# libglpk-dev + libgmp3-dev: required by igraph (clusterProfiler dep).
# python3-venv + python3-dev: for isolated scanpy/anndata env (avoids PEP 668
# "externally managed environment" block present in Ubuntu 22.04+).
# All other system libs already present in bioconductor_docker:RELEASE_3_20.
RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        libglpk-dev \
        libgmp3-dev \
        python3-venv \
        python3-dev \
    && rm -rf /var/lib/apt/lists/*

# ---- Python virtual environment + packages ------------------------------------
# Isolated venv avoids conflicts with the system Python and lets reticulate
# find a stable, fully-controlled Python via RETICULATE_PYTHON.
ENV VIRTUAL_ENV=/opt/scannex-py
RUN python3 -m venv ${VIRTUAL_ENV}
ENV PATH="${VIRTUAL_ENV}/bin:${PATH}"
ENV RETICULATE_PYTHON="${VIRTUAL_ENV}/bin/python"

RUN pip install --no-cache-dir \
    "scanpy>=1.10.0" \
    "anndata>=0.10.0" \
    "h5py>=3.8.0" \
    "numpy>=1.24.0" \
    "pandas>=2.0.0"

# ---- Bioconductor packages -----------------------------------------------------
# Installed via BiocManager to guarantee the full (non-pruned) annotation DBs.
RUN Rscript -e " \
    BiocManager::install( \
        c( \
            'clusterProfiler', \
            'enrichplot', \
            'DOSE', \
            'org.Hs.eg.db', \
            'org.Mm.eg.db', \
            'AnnotationDbi' \
        ), \
        ask    = FALSE, \
        update = FALSE \
    ) \
"

# ---- CRAN packages -------------------------------------------------------------
RUN Rscript -e " \
    install.packages( \
        c( \
            'shiny', \
            'bslib', \
            'DT', \
            'plotly', \
            'htmlwidgets', \
            'shinyWidgets', \
            'shinycssloaders', \
            'reticulate', \
            'ggplot2', \
            'ggrepel', \
            'ggridges', \
            'RColorBrewer', \
            'viridis', \
            'dplyr', \
            'tidyr', \
            'data.table', \
            'jsonlite', \
            'stringr' \
        ), \
        repos = 'https://cloud.r-project.org' \
    ) \
"

# ---- Sanity check --------------------------------------------------------------
RUN Rscript -e " \
    required <- c( \
        'shiny', 'bslib', 'DT', 'plotly', 'shinycssloaders', \
        'reticulate', 'ggplot2', 'dplyr', 'data.table', \
        'clusterProfiler', 'enrichplot', 'org.Hs.eg.db' \
    ); \
    missing <- required[!sapply(required, requireNamespace, quietly = TRUE)]; \
    if (length(missing) > 0) stop('Missing required packages: ', paste(missing, collapse=', ')); \
    message('All required packages present.') \
" \
&& python3 -c "import scanpy, anndata; print('Python: scanpy', scanpy.__version__, '| anndata', anndata.__version__)"
