# Singularity in WSL2 - Known Limitations

## Installation Status
- **Version Installed:** SingularityCE 4.1.0 (0.1.0+7-gfa01d1a-dirty)
- **Location:** `/usr/local/bin/singularity`
- **Date:** 2026-01-21

## WSL2 Limitations

### Securebits Error
When attempting to run containers in WSL2, you may encounter:
```
ERROR: Failed to set securebits: Invalid argument
```

**Root Cause:** WSL2 kernel does not fully support Linux security capabilities (securebits) that Singularity requires for privilege management.

**Impact:** 
- Singularity cannot run containers in standard mode
- SUID operations fail
- Fakeroot mode may have limited functionality

### Workarounds

#### Option 1: Use Conda Profile (RECOMMENDED for WSL2)
```bash
nextflow run main.nf -profile conda
```

**Pros:**
- Works reliably in WSL2
- No kernel limitations
- Reproducible environments

**Cons:**
- Slower first run (builds environments)
- Larger disk usage

#### Option 2: Use Docker Desktop for WSL2
If Docker Desktop is installed:
```bash
nextflow run main.nf -profile docker
```

#### Option 3: Run on Native Linux or HPC
- Singularity works best on native Linux systems
- HPC clusters typically have proper kernel support
- Production deployments should use HPC with Singularity

### Future Considerations

For production use:
1. **Development:** Use conda profile in WSL2
2. **Testing:** Test with Docker if available in WSL2
3. **Production:** Deploy to HPC cluster with native Singularity support

## Testing Strategy

Given these limitations, the testing approach is:

1. **Phase 3 (Current):** Test with conda profile in WSL2
   - Validates pipeline logic
   - Validates parameter handling
   - Validates output generation

2. **Phase 4 (Future):** Test with Singularity on HPC
   - Validates container execution
   - Validates resource management
   - Validates production readiness

## Cleanup

Build files have been removed from workspace:
- `/home/damo/scAnnex/singularity-ce-4.1.0/` (299MB) - REMOVED
- `/home/damo/scAnnex/singularity-ce-4.1.0.tar.gz` (17MB) - REMOVED

Binary remains installed at `/usr/local/bin/singularity` for future use on native Linux.
