
## Troubleshooting

### Out of Memory Errors

For large datasets (>100k cells), reduce memory usage:

```bash
nextflow run main.nf \
  --input large_data.h5ad \
  --n_top_genes 1000 \      # Reduce HVGs
  --n_pcs 30 \              # Reduce PCs
  --harmony_max_iter 5 \    # Reduce Harmony iterations
  -profile docker
```

### Slow Integration

Harmony integration can be slow on large datasets. Use optimized parameters:

```bash
nextflow run main.nf \
  --input data.h5ad \
  --run_integration true \
  --n_top_genes 1000 \
  --n_pcs 20 \
  --harmony_max_iter 5 \
  -profile docker
```

### Docker Permission Errors (Linux)

Add user to docker group:

```bash
sudo usermod -aG docker $USER
newgrp docker
```
