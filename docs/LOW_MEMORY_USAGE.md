# Running scAnnex on Low-Memory Systems (Laptops, WSL2)

## Quick Start

For running on systems with limited RAM (~4-8GB available), use the **laptop profile**:

```bash
nextflow run main.nf \
  -profile docker,laptop \
  --input samplesheet.csv \
  --outdir results
```

## Understanding the Configuration

### Default Configuration (HPC/Server)
- **Max Memory**: 128 GB
- **Max CPUs**: 16
- **Process Memory**: 6-200 GB depending on process label

### Laptop Profile (`-profile laptop`)
- **Max Memory**: 2.5 GB
- **Max CPUs**: 2
- **Process Memory**: 1.5-2.5 GB depending on process label

### How It Works

1. **Base Config** (`conf/base.config`):
   - Defines default resources for HPC environments
   - Uses `check_max()` function to automatically cap resources at `params.max_memory` and `params.max_cpus`
   - Implements retry strategy with increased resources on failure

2. **Low Memory Config** (`conf/low_memory.config`):
   - Overrides resource limits for constrained environments
   - Sets conservative memory allocations (1.5-2.5 GB per process)
   - Enables when using `-profile laptop`

3. **Profile Combination**:
   - Profiles can be combined: `-profile docker,laptop`
   - `docker` enables Docker containerization
   - `laptop` applies low-memory resource constraints

## Resource Allocation by Process Label

| Process Label          | Default (HPC) | Laptop Profile |
|------------------------|---------------|----------------|
| `process_single`       | 6 GB          | 1.5 GB         |
| `process_low`          | 12 GB         | 2 GB           |
| `process_medium`       | 36 GB         | 2.5 GB         |
| `process_high`         | 72 GB         | 2.5 GB         |
| `process_high_memory`  | 200 GB        | 2.5 GB         |

## Troubleshooting

### Still Getting Memory Errors?

If you still encounter memory errors with the laptop profile:

1. **Reduce dataset size** for testing:
   ```python
   # In Python/scanpy
   adata = adata[:5000, :]  # Use only 5000 cells
   adata.write_h5ad('test_subset.h5ad')
   ```

2. **Manually override max_memory**:
   ```bash
   nextflow run main.nf \
     -profile docker,laptop \
     --max_memory 2.GB \
     --input samplesheet.csv
   ```

3. **Check available memory in WSL2**:
   ```bash
   free -h
   ```
   
   If WSL2 has less than 4GB, you may need to configure `.wslconfig`:
   ```ini
   # %USERPROFILE%\.wslconfig (Windows)
   [wsl2]
   memory=8GB
   ```

### Process-Specific Memory Issues

If a specific process fails, check the error message for the process name and consider:

1. Using a smaller test dataset
2. Running on a system with more resources
3. Opening an issue on GitHub with the error details

## Performance Considerations

Running single-cell analysis on low-memory systems:
- **Slower execution**: Processes may swap to disk
- **Limited dataset size**: Recommended <10K cells for laptops
- **Serial processing**: Only 1-2 processes run concurrently

For production runs with large datasets (>20K cells), use:
- HPC cluster with Slurm/PBS
- Cloud instance with 32+ GB RAM
- Default profile without laptop restrictions

## Custom Resource Limits

You can also set custom limits without using the laptop profile:

```bash
nextflow run main.nf \
  -profile docker \
  --max_memory 4.GB \
  --max_cpus 2 \
  --max_time 6.h \
  --input samplesheet.csv
```

This allows fine-tuning based on your specific system capabilities.

## Monitoring Resource Usage

Track resource consumption during execution:

```bash
# Monitor in real-time
watch -n 2 'docker stats --no-stream'

# Check Nextflow execution report (after run)
ls results/pipeline_info/execution_report.html
```

The execution report includes detailed resource usage for each process, helping optimize future runs.
