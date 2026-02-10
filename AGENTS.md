# AGENTS.md

> This file is the authoritative guide for all agentic coding agents working in the scAnnex repository. It documents project setup, installation, build/lint/test commands, dependency management, and code style conventions for maintaining high code quality and reproducibility.

---

## 1. Project Overview

**scAnnex** is an automated Nextflow DSL2 pipeline for single-cell RNA-seq analysis with integrated analytical modules and R Shiny dashboard. The project supports Conda, Docker, Wave, and Singularity environments. Results can be explored interactively via the dashboard.

---

## 2. Setup & Installation (Ubuntu/Linux Best Practices)

### Prerequisite Software (required for dev & most agent tasks)
- **Nextflow** `>=23.04.0`
- **Docker** (for container profiles, dashboards, and some tests)
- **Conda** (or Mamba; Miniforge recommended)
- **Python** (`>=3.10` pipeline, `>=3.11` dashboard)
- **R** (`>=4.3`) for dashboard
- (Optional: Singularity/Apptainer for HPC)

### Explicit Install Instructions

#### 2.1. Docker (latest stable)
```bash
sudo apt-get remove docker docker-engine docker.io containerd runc || true
sudo apt-get update
sudo apt-get install -y ca-certificates curl gnupg lsb-release
sudo install -m 0755 -d /etc/apt/keyrings
curl -fsSL https://download.docker.com/linux/ubuntu/gpg | sudo gpg --dearmor -o /etc/apt/keyrings/docker.gpg
echo \ \
  "deb [arch=$(dpkg --print-architecture) signed-by=/etc/apt/keyrings/docker.gpg] https://download.docker.com/linux/ubuntu \ \
  $(lsb_release -cs) stable" | \
  sudo tee /etc/apt/sources.list.d/docker.list > /dev/null
sudo apt-get update
sudo apt-get install -y docker-ce docker-ce-cli containerd.io docker-buildx-plugin docker-compose-plugin
sudo docker run hello-world # optionally test
```

#### 2.2. Nextflow (latest stable)
```bash
curl -s https://get.nextflow.io | bash
sudo mv nextflow /usr/local/bin/
```

#### 2.3. Conda (Miniforge recommended)
```bash
curl -L -O https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
bash Miniforge3-Linux-x86_64.sh
```

#### 2.4. R and Dashboard Deps
The dashboard setup script will install all required R and Python dependencies using conda:
```bash
cd dashboard
./setup_dashboard.sh
```

#### 2.5. Python linting tools
Required for CI and code formatting:
```bash
pip install black flake8
```

---

## 3. Dependency Management
- **Pipeline** dependencies: handled by Nextflow profile (Conda or Docker).
- **Dashboard** dependencies: auto-created by `dashboard/environment_dashboard.yml` and installed by setup script.
- **Python/R deps**: Managed via conda env or containers, not manually.

---

## 4. Build, Lint & Test Commands

### 4.1. Run Pipeline
- **Demo:**
  ```bash
  nextflow run main.nf --input data_demo/10xMTX/samplesheet.csv --outdir results -profile conda
  ```
- **Custom data:**
  ```bash
  nextflow run main.nf -profile conda --input samplesheet.csv --outdir results --max_memory '8.GB'
  ```

### 4.2. Run Single Test (Pipeline module)
- Analytical core test:
  ```bash
  bash tests/test_analytical_core.sh
  ```
- Integration quick test:
  ```bash
  bash tests/test_integration_quick.sh
  ```

### 4.3. Dashboard Testing
- Run all dashboard tests (from dashboard dir):
  ```bash
  Rscript test_dashboard_full.R
  # Or partial:
  Rscript test_dashboard_data.R
  ```
- Launch dashboard interactively:
  ```bash
  ./launch_dashboard.sh ../results
  # or
  R -e "shiny::runApp('.', host='127.0.0.1', port=8888)"
  ```

### 4.4. Linting & Formatting
- Python:
  ```bash
  black --check bin/*.py
  flake8 bin/*.py --max-line-length=120 --extend-ignore=E203,W503
  ```
- Please run these before PR or as required by CI.

### 4.5. CI Workflow
- Pushes and PRs to `main` or `dev` trigger: lint, unit/integration test, Docker profile test.
- Test data: `data_demo/`.

---

## 5. Code Style Guidelines

### 5.1. Python (pipeline & scripts)
- **Formatting:** Black (`line length 88`, CI allows 120).
- **Linting:** Flake8 (`--max-line-length=120 --extend-ignore=E203,W503`).
- **Imports:** Stdlib, then external, then local (each group separated by a blank line).
- **Typing:** Prefer explicit type hints on all function signatures.
- **Naming:**
  - Functions/vars: `snake_case`
  - Classes: `CamelCase`
  - Constants: `UPPER_CASE`
- **Error handling:** Use `try`/`except` for recoverable errors; *always* log actionable errors; fail loudly for configuration/env issues.
- **Other conventions:**
  - Avoid very large functions (prefer decomposition).
  - Public functions/classes require docstrings (simple=one liner, complex=Google/numpy style).
  - Avoid unused variables/imports (CI will fail).
  - CLI entrypoints must use `if __name__ == "__main__":` guard.

### 5.2. R (dashboard)
- **Indent:** 2 spaces per level.
- **Assignment:** Use `<-` for assignment except in args.
- **Max line:** ~100 chars (preferred).
- **Naming:** `snake_case` for vars/funs, `UPPER_CASE` for constants.
- **Error handling:** Use `stop()`/`tryCatch()` as appropriate, always inform via UI if possible.
- **Structure:** Global variables/functions in `global.R`. UI in `ui.R`, server logic in `server.R`.
- **Comments:** Use `#`, doc all major/exported functions.
- **Reproducibility:** Only use packages in `environment_dashboard.yml`.

### 5.3. YAML/Nextflow
- Two-space indentation.
- Lowercase profile/environment names.
- Document non-obvious params.

---

## 6. Data Conventions (For Agents)
- **H5AD files:**
  - Must have `obs['predicted_labels']`, `obs['celltypist_score']`, `obsm['X_umap']`.
  - See `dashboard/README.md` for more.
- **Outputs:**
  - Follow strict directory/filename conventions for compatibility with downstream tools.

---

## 7. General Notes
- All development/testing should be done in an isolated conda environment or within a Docker container when possible.
- If adding new features, update both this file and README/data conventions so agents and humans remain in sync.
- There are no Cursor or Copilot rule files in this repository. If added, document and integrate them here immediately.
- Keep this file up to date when style or workflow requirements change.

---

**Maintain this file as a living resource. Agentic coding agents must always refer to, follow, and update this file in sync with repository best practices.**
