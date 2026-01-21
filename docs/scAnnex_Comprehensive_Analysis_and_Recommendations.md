# scAnnex Pipeline: Comprehensive Analysis & Recommendations

**Analysis Date:** January 21, 2026  
**Repository:** https://github.com/damouzo/scAnnex  
**Current Version:** SLC v1.0 (Simple, Lovable, Complete)  
**Analyst:** Seqera AI - Nextflow Expert System

---

## 📊 EXECUTIVE SUMMARY

### Overall Assessment: **PROMISING BUT NEEDS CRITICAL FIXES** ⚠️

**Viability Score:** 6.5/10  
**Time to Production-Ready:** 2-4 weeks (with focused effort)  
**Key Verdict:** The pipeline has solid architecture and good design philosophy, but contains **critical syntax errors** that will prevent execution. The modular structure is excellent, documentation is comprehensive, but the code needs immediate linting, testing, and compliance with modern Nextflow standards.

---

## ✅ STRENGTHS - What's Working Well

### 1. **Excellent Architecture & Design Philosophy**
- **SLC Approach (Simple, Lovable, Complete):** Smart pivot from MVP to end-to-end functionality
- **Modular Structure:** Clean separation of concerns with well-defined modules
- **Comprehensive Documentation:** TODO.md, SUMMARY.md, and multiple guides show attention to detail
- **User-Centric Features:** Cell Attrition Log is an outstanding transparency feature
- **Multi-Resolution Clustering:** Avoids single-resolution bias—this is best practice

### 2. **Strong Bioinformatics Foundation**
- **Scanpy-Based Workflow:** Industry-standard single-cell analysis
- **CellTypist Integration:** Modern, accurate auto-annotation
- **Quantile-Based QC:** More robust than fixed thresholds
- **Harmony Integration:** Gold-standard batch correction
- **Dashboard-Ready Outputs:** CSV exports (umap_coordinates.csv, cell_metadata.csv) show foresight

### 3. **Good Container Strategy**
- Docker and Singularity support
- Extended scanpy container with all dependencies
- Proper environment isolation (PYTHONNOUSERSITE, R_PROFILE_USER)

### 4. **Thoughtful Resource Management**
- Low memory profile for 8GB laptops
- Process-specific resource labels
- Configurable max_memory, max_cpus, max_time

---

## 🚨 CRITICAL ISSUES - Must Fix Before First Use

### 1. **BLOCKER: Nextflow Syntax Error in base.config**

**Location:** `conf/base.config:61`  
**Error:** `Unexpected input: '('` in function definition

**Issue:**
```groovy
// ❌ INCORRECT - Functions cannot be defined in config files (strict mode)
def check_max(obj, type) {
    // function body
}
```

**Impact:** **PIPELINE WILL NOT RUN** - This is a compilation error

**Fix Required:**
```groovy
// ✅ OPTION 1: Move to lib/Utils.groovy (RECOMMENDED)
// Create: lib/Utils.groovy
class Utils {
    static def checkMax(obj, type, params) {
        if (type == 'memory') {
            try {
                if (obj.compareTo(params.max_memory as nextflow.util.MemoryUnit) == 1)
                    return params.max_memory as nextflow.util.MemoryUnit
                else
                    return obj
            } catch (all) {
                println "   ### ERROR ###   Max memory '${params.max_memory}' is not valid! Using default value: $obj"
                return obj
            }
        } else if (type == 'time') {
            try {
                if (obj.compareTo(params.max_time as nextflow.util.Duration) == 1)
                    return params.max_time as nextflow.util.Duration
                else
                    return obj
            } catch (all) {
                println "   ### ERROR ###   Max time '${params.max_time}' is not valid! Using default value: $obj"
                return obj
            }
        } else if (type == 'cpus') {
            try {
                return Math.min(obj, params.max_cpus as int)
            } catch (all) {
                println "   ### ERROR ###   Max cpus '${params.max_cpus}' is not valid! Using default value: $obj"
                return obj
            }
        }
    }
}

// ✅ OPTION 2: Remove function and inline checks (SIMPLER)
// In conf/base.config, replace all check_max() calls:
process {
    cpus   = { Math.min(1, params.max_cpus as int) }
    memory = { 
        def mem = 6.GB * task.attempt
        mem.compareTo(params.max_memory as nextflow.util.MemoryUnit) == 1 ? 
            params.max_memory as nextflow.util.MemoryUnit : mem
    }
    // ... etc
}
```

**Priority:** 🔴 **CRITICAL - FIX IMMEDIATELY**

---

### 2. **BLOCKER: Missing Validation Plugin Configuration**

**Location:** `main.nf:14`  
**Issue:**
```groovy
include { validateParameters; paramsHelp } from 'plugin/nf-validation'
```

**Problem:** `nf-validation` plugin not declared in config or plugins section

**Impact:** Pipeline will fail with "Module not found" error

**Fix Required:**
```groovy
// Add to nextflow.config
plugins {
    id 'nf-validation@1.1.3'
}

// OR remove validation if not needed yet:
// Comment out in main.nf:
// include { validateParameters; paramsHelp } from 'plugin/nf-validation'
// if (params.help) {
//     log.info paramsHelp("nextflow run main.nf --input samplesheet.csv")
//     System.exit(0)
// }
// validateParameters()
```

**Priority:** 🔴 **CRITICAL**

---

### 3. **HIGH RISK: Input Channel Parsing Expects CSV Samplesheet**

**Location:** `workflows/scannex.nf:25-35`  
**Issue:**
```groovy
Channel
    .fromPath(params.input)
    .splitCsv(header: true)
    .map { row ->
        def meta = [
            id: row.sample_id,
            file_type: row.file_type,
            batch: row.batch,
            condition: row.condition
        ]
        [ meta, file(row.file_path, checkIfExists: true) ]
    }
```

**Problems:**
1. README.md examples show `--input sample_data.h5ad` (single file)
2. Workflow expects `--input samplesheet.csv` (CSV with multiple samples)
3. **MISMATCH between documentation and implementation**

**Impact:** All documented examples will fail

**Fix Required:**
```groovy
// ✅ Add flexible input handling
def input_ch
if (params.input.endsWith('.csv') || params.input.endsWith('.tsv')) {
    // CSV/TSV samplesheet
    input_ch = Channel
        .fromPath(params.input)
        .splitCsv(header: true, sep: params.input.endsWith('.tsv') ? '\t' : ',')
        .map { row ->
            def meta = [
                id: row.sample_id ?: file(row.file_path).simpleName,
                file_type: row.file_type ?: params.input_type,
                batch: row.batch ?: 'batch1',
                condition: row.condition ?: 'default'
            ]
            [ meta, file(row.file_path, checkIfExists: true) ]
        }
} else {
    // Single file input
    input_ch = Channel
        .fromPath(params.input, checkIfExists: true)
        .map { file ->
            def meta = [
                id: file.simpleName,
                file_type: params.input_type,
                batch: 'batch1',
                condition: 'default'
            ]
            [ meta, file ]
        }
}
```

**Priority:** 🔴 **HIGH - Documentation vs Implementation Conflict**

---

### 4. **Unused Variables Warning**

**Location:** `subworkflows/local/utils_nfcore_scannex_pipeline.nf`  
**Issue:** Multiple declared but unused variables (lines 28-31, 101-108)

**Impact:** Code clutter, potential confusion, linting warnings

**Fix:** Remove unused variables or implement functionality

**Priority:** 🟡 **MEDIUM - Code Quality**

---

### 5. **Conda Environment Path in Process Directives**

**Location:** All module files in `modules/local/*.nf`  
**Warning:** `The use of projectDir in a process is discouraged`

**Issue:**
```groovy
conda "${projectDir}/env/scanpy.yml"
```

**Problem:**
- `${projectDir}` in processes is anti-pattern in Nextflow DSL2
- Conda environment files not included in repository (env/ directory missing)

**Fix Required:**
```groovy
// ✅ OPTION 1: Use full conda specification (RECOMMENDED)
conda 'bioconda::scanpy=1.10.0 bioconda::anndata=0.10.3 conda-forge::numpy=1.24.0'

// ✅ OPTION 2: Create env/ directory and add scanpy.yml
mkdir env
# Create env/scanpy.yml with all dependencies

// ✅ OPTION 3: Remove conda directive and rely on containers
// Just delete the conda line - containers are already defined
```

**Priority:** 🟡 **MEDIUM - Will fail in conda profile**

---

## ⚠️ HIGH-PRIORITY IMPROVEMENTS

### 6. **No End-to-End Testing Evidence**

**Issue:** TODO.md shows "END-TO-END TEST" as completed ✅, but:
- No test results in repository
- No CI/CD pipeline (GitHub Actions empty)
- Test profile expects files that may not exist (`test_data/outputs/PBMC_MTX_quick_test.h5ad`)

**Risk:** Pipeline may fail on first real execution

**Recommendation:**
```yaml
# Add .github/workflows/ci.yml
name: Pipeline Tests
on: [push, pull_request]
jobs:
  test:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v3
      - name: Install Nextflow
        run: |
          wget -qO- https://get.nextflow.io | bash
          sudo mv nextflow /usr/local/bin/
      - name: Download Test Data
        run: python bin/download_test_data.py --output-dir test_data
      - name: Run Pipeline Test
        run: nextflow run main.nf -profile test,docker --outdir test_results
      - name: Verify Outputs
        run: |
          test -f test_results/qc_results/qc_report.json
          test -f test_results/standard_processing_results/umap_coordinates.csv
```

**Priority:** 🟠 **HIGH - Trust & Reliability**

---

### 7. **Dashboard Implementation Incomplete**

**Status:** Core files exist but marked "🚧 IN PROGRESS" in TODO.md

**Issues:**
- No evidence of successful dashboard launch
- No screenshots or demo
- Multiple launch scripts suggest troubleshooting issues

**Recommendation:**
1. Complete minimal viable dashboard
2. Test with real pipeline outputs
3. Record screencast demonstration
4. Add screenshots to README.md

**Priority:** 🟠 **HIGH - User Experience**

---

### 8. **Container Image Tags Not Pinned**

**Issue:**
```groovy
container "quay.io/biocontainers/scanpy:1.7.2--pyhdfd78af_0"
```

**Problems:**
- Scanpy 1.7.2 is from 2021 (outdated)
- No Wave integration despite documentation mentioning it
- Manual container management is error-prone

**Recommendation:**
```groovy
// ✅ Use Wave for automatic container generation
process QUALITY_CONTROL {
    conda 'bioconda::scanpy=1.10.0'
    // Wave will auto-generate container from conda spec
    
    // OR use your custom extended container:
    container "damouzo/scanpy-extended:1.10.0"
}
```

**Priority:** 🟠 **HIGH - Reproducibility**

---

### 9. **Missing Essential Scripts**

**Expected in `bin/` but not verified:**
- `quality_control.py` - ✅ EXISTS
- `doublet_detection.py` - ✅ EXISTS
- `standard_processing.py` - ✅ EXISTS
- `normalize_integrate.py` - ✅ EXISTS
- `auto_annot_celltypist.py` - ✅ EXISTS
- `unify_input.py` - ✅ EXISTS

**Good news:** All Python scripts exist!

**BUT:** No evidence of Python linting or testing

**Recommendation:**
```bash
# Add to repository root
make test-python:
    cd bin && python -m pytest tests/
    cd bin && flake8 *.py
    cd bin && black --check *.py
```

**Priority:** 🟡 **MEDIUM - Code Quality**

---

### 10. **Parameter Validation Schema Incomplete**

**File:** `nextflow_schema.json`  
**Issue:** Schema exists but many parameters lack:
- Descriptions
- Type validation
- Enum constraints
- Default value documentation

**Recommendation:**
```json
{
  "clustering_resolutions": {
    "type": "string",
    "default": "0.1,0.3,0.5,0.7,0.9",
    "pattern": "^[0-9]+(\\.[0-9]+)?(,[0-9]+(\\.[0-9]+)?)*$",
    "description": "Comma-separated list of Leiden clustering resolutions",
    "help_text": "Higher values = more fine-grained clusters"
  }
}
```

**Priority:** 🟡 **MEDIUM - User Experience**

---

## 🔮 OBSTACLES & RISKS

### 1. **First-Time User Experience Will Be Frustrating**

**Why:**
- Syntax errors will prevent execution
- Documentation examples don't match implementation
- No verified working example in README

**Impact:** High bounce rate, negative reputation

**Mitigation:**
- Fix all CRITICAL issues before public announcement
- Add "Quick Start in 30 Seconds" with verified commands
- Record video tutorial

---

### 2. **Conda Profile is Broken**

**Why:**
- `env/scanpy.yml` file doesn't exist
- `projectDir` usage in process directives
- No conda testing evidence

**Impact:** Users preferring conda over Docker will fail

**Mitigation:**
- Create complete conda environment file
- Test conda profile separately
- Document limitations (recommend Docker first)

---

### 3. **Resource Requirements Unclear**

**Current Documentation:**
- Says "works on 8GB laptops" but test profile sets `max_memory = '8.GB'`
- No benchmarking data (e.g., "1k cells = 2GB RAM, 5 minutes")
- No guidance on when to use HPC vs local

**Impact:** Users waste time on undersized systems

**Mitigation:**
```markdown
## Hardware Requirements

| Dataset Size | Recommended RAM | Estimated Time | Profile |
|--------------|----------------|----------------|---------|
| < 5,000 cells | 8 GB | 10-15 min | `docker,laptop` |
| 5,000-50,000 | 16 GB | 30-60 min | `docker` |
| 50,000-200,000 | 32 GB | 2-4 hours | `singularity` (HPC) |
| > 200,000 | 64+ GB | 4-12 hours | `singularity` (HPC) |
```

---

### 4. **No Community Yet**

**Risks:**
- No GitHub stars/forks suggests no external testing
- No issues/discussions = no feedback loop
- Solo development = knowledge silo

**Mitigation:**
- Announce on Seqera Community Slack
- Post to nf-core forum
- Create Twitter/LinkedIn posts with demo

---

### 5. **Dashboard Complexity May Delay Release**

**Risk:** Dashboard has 28 files, multiple launch methods, troubleshooting docs

**This suggests:**
- Significant technical challenges
- Platform-specific issues (WSL2 troubleshooting)
- Not yet production-ready

**Mitigation:**
- Release pipeline WITHOUT dashboard first (v0.9)
- Mark dashboard as "beta" feature
- Focus on CSV outputs being usable in RStudio/Python notebooks

---

## 📅 REALISTIC TIMELINE TO PRODUCTION

### Week 1: **Critical Fixes** (5-8 hours)
- [ ] Fix `check_max()` function in base.config
- [ ] Add nf-validation plugin OR remove validation calls
- [ ] Fix input channel to accept single files AND samplesheets
- [ ] Test with `nextflow lint` and fix all errors
- [ ] Run full pipeline on test data and document success

### Week 2: **Validation & Testing** (10-15 hours)
- [ ] Set up GitHub Actions CI/CD
- [ ] Run pipeline on 3 different datasets (small, medium, large)
- [ ] Create benchmark table (cells → time/memory)
- [ ] Fix all bugs discovered during testing
- [ ] Create conda environment files or remove conda profile

### Week 3: **Documentation & Polish** (8-12 hours)
- [ ] Update README with working examples (verified)
- [ ] Create quickstart video (5 minutes max)
- [ ] Add troubleshooting section with real issues
- [ ] Pin all container versions
- [ ] Add CITATION.cff file

### Week 4: **Dashboard & Release** (15-20 hours)
- [ ] Complete dashboard testing
- [ ] Create dashboard demo video
- [ ] Write release notes
- [ ] Announce on community platforms
- [ ] Monitor first-week user issues

**Total Estimated Effort:** 38-55 hours  
**Calendar Time:** 2-4 weeks (depending on availability)

---

## 🎯 IMMEDIATE ACTION ITEMS (Next 24-48 Hours)

### Priority 1: Make It Run 🔴

```bash
# 1. Fix base.config
# Create lib/Utils.groovy and move check_max() function
# OR inline all checks

# 2. Fix main.nf
# Remove or fix nf-validation plugin dependency

# 3. Fix input handling
# Accept both single files and samplesheets

# 4. Test
nextflow run main.nf --help  # Should not error
nextflow lint .              # Should have 0 errors
```

### Priority 2: Verify Test Works 🟠

```bash
# 5. Run test
python bin/download_test_data.py --output-dir test_data
nextflow run main.nf \
  --input test_data/outputs/PBMC_MTX_quick_test.h5ad \
  --input_type h5ad \
  --outdir test_results \
  -profile docker

# 6. Verify outputs exist
ls test_results/qc_results/
ls test_results/standard_processing_results/
```

### Priority 3: Document Success 🟡

```bash
# 7. Update README with VERIFIED commands
# 8. Take screenshots of outputs
# 9. Commit and push
git add -A
git commit -m "Fix critical syntax errors, verified working test"
git push
```

---

## 📝 DETAILED FIX CHECKLIST

### For OpenCode or Development Team

#### 🔴 CRITICAL - Do These First

- [ ] **FIX conf/base.config:61**
  - Option A: Create `lib/Utils.groovy` and move `check_max()` function
  - Option B: Inline all `check_max()` calls in base.config
  - Verify with `nextflow lint conf/base.config`

- [ ] **FIX main.nf validation plugin**
  - Option A: Add `plugins { id 'nf-validation' }` to nextflow.config
  - Option B: Remove validation calls from main.nf
  - Test with `nextflow run main.nf --help`

- [ ] **FIX workflows/scannex.nf input channel**
  - Add logic to detect single file vs CSV samplesheet
  - Support both input methods
  - Update README examples to match

- [ ] **RUN nextflow lint and fix ALL errors**
  ```bash
  nextflow lint .
  # Fix every error until output shows "✅ All files had no errors"
  ```

- [ ] **TEST with actual data**
  ```bash
  python bin/download_test_data.py --output-dir test_data
  nextflow run main.nf --input <path> --outdir test_results -profile docker
  # Must complete without errors
  ```

#### 🟠 HIGH Priority - Do This Week

- [ ] **Create env/scanpy.yml for conda profile**
  ```yaml
  name: scannex-scanpy
  channels:
    - conda-forge
    - bioconda
  dependencies:
    - python=3.10
    - scanpy=1.10.0
    - anndata=0.10.3
    - celltypist=1.6.2
    - scrublet=0.2.3
    - harmonypy=0.0.9
    # ... all dependencies
  ```

- [ ] **Pin container versions**
  - Replace `scanpy:1.7.2` with `scanpy:1.10.0`
  - Document container source
  - Consider building and pushing to Docker Hub

- [ ] **Add GitHub Actions CI**
  ```yaml
  .github/workflows/ci.yml:
    - Test pipeline with small dataset
    - Verify outputs exist
    - Fail if any errors occur
  ```

- [ ] **Fix unused variables in subworkflows/local/utils_nfcore_scannex_pipeline.nf**
  - Remove or implement unused variables

- [ ] **Complete dashboard testing**
  - Launch dashboard with real outputs
  - Document any errors and fix
  - Take screenshots for README

#### 🟡 MEDIUM Priority - This Month

- [ ] **Improve nextflow_schema.json**
  - Add descriptions for all parameters
  - Add validation patterns
  - Add help_text with examples

- [ ] **Add Python testing**
  ```bash
  cd bin
  pytest tests/
  flake8 *.py
  black --check *.py
  ```

- [ ] **Create benchmarking table**
  - Test 1k, 10k, 100k cell datasets
  - Record time and memory usage
  - Add to README

- [ ] **Verify all documentation links work**
  - Check all docs/*.md files
  - Fix broken cross-references
  - Remove outdated information

- [ ] **Create CITATION.cff**
  ```yaml
  cff-version: 1.2.0
  title: scAnnex
  message: "If you use this software, please cite it as below."
  authors:
    - family-names: "Your Name"
  version: 1.0.0
  date-released: 2026-XX-XX
  ```

#### 🟢 LOW Priority - When Time Permits

- [ ] **Add nf-core template compliance**
  - Use `nf-core create` to generate structure
  - Adopt nf-core modules where applicable

- [ ] **Implement Wave integration**
  ```groovy
  wave {
      enabled = true
      strategy = ['conda']
  }
  ```

- [ ] **Add MultiQC report**
  - Aggregate QC metrics
  - Generate HTML report

- [ ] **Create tutorial dataset**
  - Small, fast-running example
  - Include in repository or release

- [ ] **Set up Seqera Platform Cloud**
  - Create public workspace
  - Add example runs
  - Enable for community testing

---

## 💡 RECOMMENDATIONS FOR NEXT STEPS

### Immediate (This Week)

1. **Fix all CRITICAL issues** - pipeline must compile and run
2. **Test end-to-end** - verify with real data
3. **Update documentation** - make examples match implementation
4. **Enable GitHub Actions** - automate testing

### Short-term (This Month)

5. **Complete dashboard** - make it production-ready or mark as beta
6. **Add benchmarks** - give users realistic expectations
7. **Pin dependencies** - ensure reproducibility
8. **Gather feedback** - share with 3-5 beta users

### Long-term (Next Quarter)

9. **Publish paper/preprint** - describe methodology
10. **Submit to nf-core** - join the community
11. **Add advanced features** - trajectory analysis, cell-cell communication
12. **Scale testing** - validate on >100k cell datasets

---

## 🎬 CONCLUSION

### The Good News ✅

scAnnex has **excellent bones**:
- Smart architecture (SLC philosophy)
- Thoughtful features (Cell Attrition Log)
- Comprehensive documentation
- Modern tools (Scanpy, CellTypist, Harmony)

### The Bad News ❌

It's currently **broken**:
- Won't compile due to syntax errors
- Documentation doesn't match implementation
- No verified end-to-end test
- Dashboard incomplete

### The Verdict 🎯

**With 2-4 weeks of focused work**, this can become a **highly valuable community resource**. The foundation is solid; it just needs debugging, testing, and polish.

### Is It Worth Using Today? 🤔

**NO** - Not yet. Wait for:
1. Critical fixes to be merged
2. GitHub Actions showing passing tests
3. At least one verified successful run documented

### Should You Clone and Try? 🔬

**YES** - If you're willing to:
- Fix the syntax errors yourself
- Contribute fixes back to the project
- Help test and provide feedback

### Recommended Path Forward 🛤️

1. **Authors:** Fix CRITICAL issues this week
2. **Community:** Watch/star the repo, wait for v1.0 release
3. **Brave testers:** Clone, fix, test, submit PRs
4. **Everyone else:** Wait 2-4 weeks for stable release

---

## 📞 SUPPORT FOR AUTHORS

If you need help fixing these issues:

1. **Seqera Community Slack:** Post in #nextflow-help
2. **nf-core community:** Great for best practices
3. **Seqera AI:** Can help debug specific errors
4. **GitHub Discussions:** Enable for community Q&A

**This pipeline has real potential. Don't let syntax errors hold it back!**

---

**Analysis Completed:** January 21, 2026  
**Analyst:** Seqera AI - Nextflow DSL2 Expert  
**Report Version:** 1.0  
**Next Review:** After critical fixes implemented

---

## 🔗 Quick Reference Links

- **Repository:** https://github.com/damouzo/scAnnex
- **Nextflow Linting Docs:** https://nextflow.io/docs/latest/cli.html#lint
- **nf-validation Plugin:** https://nextflow-io.github.io/nf-validation/
- **Wave Containers:** https://seqera.io/wave/
- **Seqera Community:** https://community.seqera.io

---

**End of Analysis Report**
