#!/usr/bin/env python3
"""
Download and prepare a small public multi-sample demo dataset.

Source dataset:
- GEO GSE96583 (Kang et al., IFN-beta PBMC)
- Uses batch 2 control (GSM2560248) and stim (GSM2560249)

Format-specific outputs are created under data_demo/MultiSample/<FORMAT>/:
- H5AD: 4 H5AD files + samplesheet.csv + contrasts_example.csv
- MTX: 4 10x-like directories + samplesheet.csv + contrasts_example.csv
- RDS: 4 RDS files + samplesheet.csv + contrasts_example.csv
"""

import argparse
import gzip
import logging
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path
from urllib.request import urlretrieve

import anndata as ad
import numpy as np
import pandas as pd
from scipy.io import mmread, mmwrite
from scipy.sparse import csr_matrix

if hasattr(ad, "settings") and hasattr(ad.settings, "allow_write_nullable_strings"):
    ad.settings.allow_write_nullable_strings = True

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
logger = logging.getLogger(__name__)

SERIES_BASE_URL = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE96nnn/GSE96583/suppl"
GSM_CONTROL_BASE_URL = "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2560nnn/GSM2560248/suppl"
GSM_TREATED_BASE_URL = "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM2560nnn/GSM2560249/suppl"

DOWNLOAD_URLS = {
    "control_mtx": f"{GSM_CONTROL_BASE_URL}/GSM2560248_2.1.mtx.gz",
    "control_barcodes": f"{GSM_CONTROL_BASE_URL}/GSM2560248_barcodes.tsv.gz",
    "treated_mtx": f"{GSM_TREATED_BASE_URL}/GSM2560249_2.2.mtx.gz",
    "treated_barcodes": f"{GSM_TREATED_BASE_URL}/GSM2560249_barcodes.tsv.gz",
    "genes": f"{SERIES_BASE_URL}/GSE96583_batch2.genes.tsv.gz",
}


def download_file(url: str, destination: Path) -> None:
    """Download URL to destination if missing."""
    if destination.exists():
        logger.info("File already exists, skipping download: %s", destination)
        return

    destination.parent.mkdir(parents=True, exist_ok=True)
    logger.info("Downloading %s", url)
    try:
        urlretrieve(url, destination)
    except Exception as exc:
        logger.error("Failed to download %s: %s", url, exc)
        sys.exit(1)


def read_genes(genes_file: Path) -> pd.DataFrame:
    """Read gene table from GEO supplementary TSV."""
    genes = pd.read_csv(genes_file, sep="\t", header=None)
    if genes.shape[1] == 1:
        genes.columns = ["gene_id"]
        genes["gene_name"] = genes["gene_id"]
    else:
        genes = genes.iloc[:, :2]
        genes.columns = ["gene_id", "gene_name"]
    return genes


def read_barcodes(barcodes_file: Path) -> pd.Series:
    """Read barcode list."""
    barcodes = pd.read_csv(barcodes_file, sep="\t", header=None).iloc[:, 0]
    return barcodes.astype(str)


def load_condition_adata(
    mtx_file: Path,
    barcodes_file: Path,
    genes_file: Path,
    sample_prefix: str,
    condition: str,
) -> ad.AnnData:
    """Load one condition matrix into AnnData."""
    logger.info("Loading matrix: %s", mtx_file.name)
    with gzip.open(mtx_file, "rb") as handle:
        matrix = mmread(handle)

    matrix = csr_matrix(matrix)
    genes = read_genes(genes_file)
    barcodes = read_barcodes(barcodes_file)

    # GEO matrix is genes x cells, convert to cells x genes for AnnData.
    if matrix.shape[0] == len(genes) and matrix.shape[1] == len(barcodes):
        matrix = matrix.T.tocsr()
    elif matrix.shape[0] == len(barcodes) and matrix.shape[1] == len(genes):
        matrix = matrix.tocsr()
    else:
        raise ValueError(
            f"Unexpected matrix dimensions {matrix.shape}; genes={len(genes)}, barcodes={len(barcodes)}"
        )

    obs = pd.DataFrame(index=[f"{sample_prefix}_{bc}" for bc in barcodes.tolist()])
    obs["condition"] = condition

    var = pd.DataFrame(index=genes["gene_name"].astype(str).tolist())
    var["gene_id"] = genes["gene_id"].astype(str).tolist()

    adata = ad.AnnData(X=matrix, obs=obs, var=var)
    adata.var_names_make_unique()
    return adata


def balanced_downsample(adata: ad.AnnData, n_cells: int, seed: int) -> ad.AnnData:
    """Subsample cells without replacement."""
    if adata.n_obs <= n_cells:
        return adata.copy()
    rng = np.random.default_rng(seed)
    idx = rng.choice(adata.n_obs, size=n_cells, replace=False)
    return adata[idx, :].copy()


def split_into_two_samples(
    adata: ad.AnnData,
    sample_a: str,
    sample_b: str,
    condition: str,
    seed: int,
) -> tuple[ad.AnnData, ad.AnnData]:
    """Split one condition AnnData into two pseudo-replicates."""
    rng = np.random.default_rng(seed)
    indices = rng.permutation(adata.n_obs)
    cut = adata.n_obs // 2
    idx_a = indices[:cut]
    idx_b = indices[cut:]

    adata_a = adata[idx_a, :].copy()
    adata_b = adata[idx_b, :].copy()

    adata_a.obs["sample_id"] = sample_a
    adata_a.obs["batch"] = "batch1"
    adata_a.obs["condition"] = condition

    adata_b.obs["sample_id"] = sample_b
    adata_b.obs["batch"] = "batch2"
    adata_b.obs["condition"] = condition

    return adata_a, adata_b


def write_contrasts(outdir: Path) -> None:
    """Write a basic contrasts file in target output directory."""
    contrasts = pd.DataFrame(
        [
            {
                "contrast_id": "treated_vs_control",
                "variable": "condition",
                "group1": "treated",
                "group2": "control",
            }
        ]
    )
    contrasts.to_csv(outdir / "contrasts_example.csv", index=False)


def write_samplesheet(outdir: Path, file_type: str, path_map: dict[str, str]) -> None:
    """Write samplesheet.csv for one output format."""
    samplesheet = pd.DataFrame(
        [
            {
                "sample_id": "control_batch1",
                "file_type": file_type,
                "file_path": path_map["control_batch1"],
                "batch": "batch1",
                "condition": "control",
            },
            {
                "sample_id": "control_batch2",
                "file_type": file_type,
                "file_path": path_map["control_batch2"],
                "batch": "batch2",
                "condition": "control",
            },
            {
                "sample_id": "treated_batch1",
                "file_type": file_type,
                "file_path": path_map["treated_batch1"],
                "batch": "batch1",
                "condition": "treated",
            },
            {
                "sample_id": "treated_batch2",
                "file_type": file_type,
                "file_path": path_map["treated_batch2"],
                "batch": "batch2",
                "condition": "treated",
            },
        ]
    )
    samplesheet.to_csv(outdir / "samplesheet.csv", index=False)


def write_h5ad_outputs(output_root: Path, samples: dict[str, ad.AnnData]) -> None:
    """Write H5AD sample files and metadata files."""
    outdir = output_root / "H5AD"
    outdir.mkdir(parents=True, exist_ok=True)

    path_map: dict[str, str] = {}
    for sample_id, sample_adata in samples.items():
        output_file = outdir / f"{sample_id}.h5ad"
        logger.info("Writing %s cells to %s", sample_adata.n_obs, output_file)
        sample_adata.write_h5ad(output_file)
        path_map[sample_id] = str(output_file)

    write_samplesheet(outdir, "h5ad", path_map)
    write_contrasts(outdir)


def write_one_mtx_sample(sample_dir: Path, sample_adata: ad.AnnData) -> None:
    """Write one sample in 10x-style MTX layout expected by unify_input.py."""
    sample_dir.mkdir(parents=True, exist_ok=True)

    matrix = sample_adata.X.tocsc().T if hasattr(sample_adata.X, "tocsc") else csr_matrix(sample_adata.X).T
    mmwrite(sample_dir / "matrix.mtx", matrix)

    with open(sample_dir / "barcodes.tsv", "w", encoding="utf-8") as handle:
        for barcode in sample_adata.obs_names.astype(str):
            handle.write(f"{barcode}\n")

    gene_ids = (
        sample_adata.var["gene_id"].astype(str).tolist()
        if "gene_id" in sample_adata.var.columns
        else sample_adata.var_names.astype(str).tolist()
    )
    gene_names = sample_adata.var_names.astype(str).tolist()
    with open(sample_dir / "features.tsv", "w", encoding="utf-8") as handle:
        for gid, gname in zip(gene_ids, gene_names):
            handle.write(f"{gid}\t{gname}\tGene Expression\n")


def write_mtx_outputs(output_root: Path, samples: dict[str, ad.AnnData]) -> None:
    """Write MTX sample directories and metadata files."""
    outdir = output_root / "MTX"
    outdir.mkdir(parents=True, exist_ok=True)

    path_map: dict[str, str] = {}
    for sample_id, sample_adata in samples.items():
        sample_dir = outdir / sample_id
        logger.info("Writing %s cells to %s", sample_adata.n_obs, sample_dir)
        write_one_mtx_sample(sample_dir, sample_adata)
        path_map[sample_id] = str(sample_dir)

    write_samplesheet(outdir, "mtx", path_map)
    write_contrasts(outdir)


def write_rds_outputs(output_root: Path, samples: dict[str, ad.AnnData], rscript_path: str) -> None:
    """Write RDS sample files by converting temporary H5AD files with existing R script."""
    outdir = output_root / "RDS"
    outdir.mkdir(parents=True, exist_ok=True)

    converter = Path("bin/convert_h5ad_to_rds.R")
    if not converter.exists():
        raise FileNotFoundError(f"R converter not found: {converter}")

    path_map: dict[str, str] = {}
    for sample_id, sample_adata in samples.items():
        with tempfile.NamedTemporaryFile(suffix=".h5ad", delete=False) as tmp:
            tmp_h5ad = Path(tmp.name)

        try:
            sample_adata.write_h5ad(tmp_h5ad)
            rds_out = outdir / f"{sample_id}.rds"
            cmd = [
                rscript_path,
                str(converter),
                "--input",
                str(tmp_h5ad),
                "--output",
                str(rds_out),
            ]
            logger.info("Converting %s to RDS", sample_id)
            result = subprocess.run(cmd, capture_output=True, text=True, check=False)
            if result.returncode != 0:
                raise RuntimeError(
                    f"RDS conversion failed for {sample_id}: {result.stderr.strip() or result.stdout.strip()}"
                )
            path_map[sample_id] = str(rds_out)
        finally:
            if tmp_h5ad.exists():
                tmp_h5ad.unlink()

    write_samplesheet(outdir, "rds", path_map)
    write_contrasts(outdir)


def main() -> None:
    parser = argparse.ArgumentParser(description="Download and prepare a public MultiSample demo dataset")
    parser.add_argument(
        "--output-dir",
        type=str,
        default="data_demo/MultiSample",
        help="Output directory root (default: data_demo/MultiSample)",
    )
    parser.add_argument(
        "--format",
        type=str,
        choices=["h5ad", "mtx", "rds"],
        default="h5ad",
        help="Output format to generate (default: h5ad)",
    )
    parser.add_argument(
        "--max-cells-per-condition",
        type=int,
        default=1200,
        help="Maximum cells to keep per condition before split (default: 1200)",
    )
    parser.add_argument("--seed", type=int, default=42, help="Random seed")
    parser.add_argument(
        "--rscript-path",
        type=str,
        default="Rscript",
        help="Rscript executable for --format rds (default: Rscript)",
    )
    parser.add_argument(
        "--keep-downloads",
        action="store_true",
        help="Keep downloaded raw GEO files in output-dir/raw_downloads",
    )
    args = parser.parse_args()

    output_root = Path(args.output_dir)
    raw_dir = output_root / "raw_downloads"
    raw_dir.mkdir(parents=True, exist_ok=True)

    local_files = {}
    for key, url in DOWNLOAD_URLS.items():
        filename = url.rsplit("/", 1)[-1]
        destination = raw_dir / filename
        download_file(url, destination)
        local_files[key] = destination

    control = load_condition_adata(
        local_files["control_mtx"],
        local_files["control_barcodes"],
        local_files["genes"],
        sample_prefix="ctrl",
        condition="control",
    )
    treated = load_condition_adata(
        local_files["treated_mtx"],
        local_files["treated_barcodes"],
        local_files["genes"],
        sample_prefix="stim",
        condition="treated",
    )

    target_cells = min(args.max_cells_per_condition, control.n_obs, treated.n_obs)
    logger.info("Balanced subsampling to %s cells per condition", target_cells)

    control = balanced_downsample(control, target_cells, args.seed)
    treated = balanced_downsample(treated, target_cells, args.seed + 1)

    control_1, control_2 = split_into_two_samples(
        control,
        sample_a="control_batch1",
        sample_b="control_batch2",
        condition="control",
        seed=args.seed,
    )
    treated_1, treated_2 = split_into_two_samples(
        treated,
        sample_a="treated_batch1",
        sample_b="treated_batch2",
        condition="treated",
        seed=args.seed + 7,
    )

    samples = {
        "control_batch1": control_1,
        "control_batch2": control_2,
        "treated_batch1": treated_1,
        "treated_batch2": treated_2,
    }

    if args.format == "h5ad":
        write_h5ad_outputs(output_root, samples)
    elif args.format == "mtx":
        write_mtx_outputs(output_root, samples)
    elif args.format == "rds":
        write_rds_outputs(output_root, samples, args.rscript_path)

    if not args.keep_downloads and raw_dir.exists():
        shutil.rmtree(raw_dir)

    logger.info("Done. Generated %s demo in %s", args.format.upper(), output_root)


if __name__ == "__main__":
    main()
