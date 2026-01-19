# scAnnex Interactive Dashboard
# Based on rocker/shiny with Python + scanpy for reading H5AD files

FROM rocker/shiny:4.3.3

LABEL maintainer="scAnnex Development Team"
LABEL description="Interactive Shiny dashboard for scRNA-seq analysis results"

# Install system dependencies
RUN apt-get update && apt-get install -y \
    python3 \
    python3-pip \
    python3-dev \
    libhdf5-dev \
    libxml2-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libpng-dev \
    git \
    && rm -rf /var/lib/apt/lists/*

# Install Python packages for reading H5AD files
RUN pip3 install --no-cache-dir \
    numpy \
    pandas \
    scanpy \
    anndata \
    h5py

# Install R packages
RUN R -e "install.packages(c( \
    'reticulate', \
    'plotly', \
    'DT', \
    'shinydashboard', \
    'shinyWidgets', \
    'viridis', \
    'RColorBrewer', \
    'data.table', \
    'jsonlite' \
), repos='https://cloud.r-project.org/')"

# Set up Python for reticulate
RUN R -e "library(reticulate); use_python('/usr/bin/python3', required = TRUE)"

# Copy dashboard files
COPY app.R /srv/shiny-server/
COPY global.R /srv/shiny-server/
COPY server.R /srv/shiny-server/
COPY ui.R /srv/shiny-server/
COPY modules/ /srv/shiny-server/modules/
COPY www/ /srv/shiny-server/www/

# Set permissions
RUN chown -R shiny:shiny /srv/shiny-server

# Expose Shiny port
EXPOSE 3838

# Run Shiny Server
CMD ["/usr/bin/shiny-server"]
