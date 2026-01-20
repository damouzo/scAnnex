# Docker Containers for scAnnex

This directory contains Dockerfiles and build scripts for custom containers used in the scAnnex pipeline.

## Available Dockerfiles

### Dockerfile.scanpy-extended

Extends the base Scanpy biocontainer with additional single-cell tools:
- Base: `quay.io/biocontainers/scanpy:1.7.2`
- Adds: scrublet, harmonypy, celltypist

**Build:**
```bash
./build-scanpy-extended.sh
# or
docker build -t scannex/scanpy-extended:1.7.2 -f Dockerfile.scanpy-extended .
```

**Use in pipeline:**
Update module files to use `scannex/scanpy-extended:1.7.2`

## Why Custom Containers?

The official Biocontainers repository:
1. Only has Scanpy up to version 1.7.2 (not 1.9.x)
2. Doesn't provide mulled containers with scrublet+scanpy or harmony+scanpy
3. Individual tool containers may miss dependencies

Our custom container solves these issues by bundling all required tools in one image.

## Alternative: Use Conda

If you don't want to build containers, use the conda profile:
```bash
nextflow run main.nf -profile conda,laptop --input samplesheet.csv
```

This automatically installs all dependencies without Docker.
