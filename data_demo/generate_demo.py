#!/usr/bin/env python3
"""
Generate synthetic PBMC demo data for scAnnex.

Creates 6 H5AD files (3 Control + 3 Treatment) mimicking output from
nf-core/scrnaseq (CellRanger → mtx_conversions). The Treatment condition
simulates IFN-beta stimulation, producing visible DE signal in volcano plots
and GSEA enrichment in interferon response pathways.

Usage (from project root):
    python data_demo/generate_demo.py

Output:
    data_demo/Ctrl_PBMC_C{1,2,3}_filtered_matrix.h5ad
    data_demo/Treat_PBMC_T{1,2,3}_filtered_matrix.h5ad
"""

import argparse
import logging
import sys
from pathlib import Path

try:
    import anndata as ad
    import numpy as np
    import pandas as pd
    import scipy.sparse as sp
except ImportError as e:
    print(f"ERROR: {e}")
    print("Run this script inside the scAnnex conda environment:")
    print("  conda activate scannex   # or the full path to your env")
    print("  python data_demo/generate_demo.py")
    sys.exit(1)

# anndata >= 0.11 compatibility
try:
    ad.settings.allow_write_nullable_strings = True
except AttributeError:
    pass

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Gene panels
# ---------------------------------------------------------------------------

# Cell-type marker genes (used to build the expression signature)
CELL_TYPE_MARKERS = {
    "T_CD4": ["CD3D", "CD3E", "CD4", "IL7R", "TRAC", "TRBC1", "CD2", "CD27", "CCR7", "LEF1"],
    "T_CD8": ["CD3D", "CD3E", "CD8A", "CD8B", "GZMK", "TRAC", "NKG7", "EOMES", "CCL5", "GZMA"],
    "B_cell": ["MS4A1", "CD79A", "CD79B", "IGHM", "BANK1", "CD22", "FCER2", "BLK", "FCRL1", "PAX5"],
    "NK_cell": ["GNLY", "NKG7", "KLRB1", "NCAM1", "FCGR3A", "KLRD1", "GZMB", "FGFBP2", "CX3CR1", "TYROBP"],
    "Monocyte_CD14": ["LYZ", "CST3", "CD14", "S100A8", "S100A9", "FCN1", "VCAN", "CTSS", "GRN", "MNDA"],
    "Monocyte_CD16": ["FCGR3A", "LST1", "AIF1", "CDKN1C", "HES4", "MS4A7", "CX3CR1", "RHOC", "SERPINA1", "LILRB2"],
}

# IFN-alpha response genes (HALLMARK_INTERFERON_ALPHA_RESPONSE subset)
IFN_ALPHA_GENES = [
    "MX1", "MX2", "ISG15", "ISG20", "IFIT1", "IFIT2", "IFIT3", "IFIT5",
    "IFI44", "IFI44L", "IFI6", "IFI27", "IFITM1", "IFITM2", "IFITM3",
    "OAS1", "OAS2", "OAS3", "OASL", "RSAD2", "HERC5", "HERC6",
    "CXCL10", "IRF7", "EIF2AK2", "TRIM22", "TRIM5", "BST2",
]

# IFN-gamma response genes (HALLMARK_INTERFERON_GAMMA_RESPONSE subset)
IFN_GAMMA_GENES = [
    "STAT1", "TAP1", "TAP2", "TAPBP", "B2M", "IRF1", "IRF9",
    "HLA-A", "HLA-B", "HLA-C", "HLA-E", "HLA-F",
    "CXCL9", "CXCL10", "CXCL11",
    "GBP1", "GBP2", "GBP4", "GBP5",
    "PSMB8", "PSMB9", "PSME1", "PSME2",
]

# All DE genes (upregulated in Treatment)
IFN_DE_GENES = list(dict.fromkeys(IFN_ALPHA_GENES + IFN_GAMMA_GENES))  # deduplicated, ordered

# Mitochondrial protein-coding genes (chromosome MT, prefix "MT-").
# The QC module flags genes starting with "MT-" to compute pct_counts_mt.
# Using 13 protein-coding MT genes with higher expression produces a
# realistic MT% distribution (~5-20%) for informative QC violin plots.
MT_GENES = [
    "MT-ND1", "MT-ND2", "MT-ND3", "MT-ND4", "MT-ND4L", "MT-ND5", "MT-ND6",
    "MT-CO1", "MT-CO2", "MT-CO3",
    "MT-ATP6", "MT-ATP8",
    "MT-CYB",
]

# Background genes — real HGNC symbols with Entrez IDs, covering diverse KEGG/GO pathways.
# All genes here map to org.Hs.eg.db and span ribosome, metabolism, signaling, splicing, etc.,
# giving a realistic gene universe for KEGG and GO enrichment analyses.
BACKGROUND_GENES = [
    # Ribosomal proteins (real HGNC, all map to Entrez)
    "RPL3", "RPL4", "RPL5", "RPL6", "RPL7", "RPL7A", "RPL8", "RPL9", "RPL10", "RPL10A",
    "RPL11", "RPL12", "RPL13", "RPL14", "RPL15", "RPL17", "RPL18", "RPL18A", "RPL19", "RPL21",
    "RPL22", "RPL23", "RPL23A", "RPL24", "RPL26", "RPL27", "RPL27A", "RPL28", "RPL29", "RPL30",
    "RPL31", "RPL32", "RPL34", "RPL35", "RPL35A", "RPL36", "RPL37", "RPL37A", "RPL38", "RPL39",
    "RPL41", "RPLP0", "RPLP1", "RPLP2",
    "RPS2", "RPS3", "RPS3A", "RPS4X", "RPS5", "RPS6", "RPS7", "RPS8", "RPS9", "RPS10",
    "RPS11", "RPS12", "RPS13", "RPS14", "RPS15", "RPS15A", "RPS16", "RPS17", "RPS19", "RPS20",
    "RPS21", "RPS23", "RPS24", "RPS25", "RPS26", "RPS27", "RPS28", "RPS29",
    # Translation factors
    "EIF1", "EIF1A", "EIF2S1", "EIF2S2", "EIF2S3", "EIF3A", "EIF3B", "EIF3C", "EIF3D",
    "EIF3E", "EIF3F", "EIF3G", "EIF3H", "EIF3I", "EIF3J", "EIF4A1", "EIF4A2", "EIF4B",
    "EIF4E", "EIF4G1", "EIF4H", "EIF5", "EIF5A", "EIF6",
    "EEF1A1", "EEF1B2", "EEF1D", "EEF1G", "EEF2",
    # mRNA splicing / hnRNPs (all real)
    "HNRNPA1", "HNRNPA2B1", "HNRNPC", "HNRNPD", "HNRNPE", "HNRNPF", "HNRNPH1", "HNRNPH2",
    "HNRNPH3", "HNRNPK", "HNRNPL", "HNRNPM", "HNRNPR", "HNRNPU", "HNRNPAB",
    "SRSF1", "SRSF2", "SRSF3", "SRSF5", "SRSF6", "SRSF7", "SRSF9",
    "SF3A1", "SF3A2", "SF3A3", "SF3B1", "SF3B2", "SF3B3", "SF3B4",
    "SNRPA", "SNRPB", "SNRPC", "SNRPD1", "SNRPD2", "SNRPD3", "SNRPE", "SNRPF", "SNRPG",
    # Ubiquitin / proteasome
    "UBB", "UBC", "UBA52", "UBL4A",
    "PSMA1", "PSMA2", "PSMA3", "PSMA4", "PSMA5", "PSMA6", "PSMA7",
    "PSMB1", "PSMB2", "PSMB3", "PSMB4", "PSMB5", "PSMB6", "PSMB7",
    "PSMC1", "PSMC2", "PSMC3", "PSMC4", "PSMC5", "PSMC6",
    "PSMD1", "PSMD2", "PSMD3", "PSMD4", "PSMD6", "PSMD7", "PSMD8",
    # Cytoskeleton
    "ACTB", "ACTG1", "TUBA1A", "TUBA1B", "TUBA4A", "TUBB", "TUBB4B", "TUBB6",
    "VIM", "FLNA", "FLNB", "MYH9", "MYL6", "MYL12A", "MYL12B",
    "CFL1", "CFL2", "ARPC1B", "ARPC2", "ARPC3", "ARPC4", "ARPC5",
    # Metabolism (glycolysis / OXPHOS)
    "GAPDH", "PGAM1", "PKM", "LDHA", "LDHB", "ENO1", "TPI1", "PGK1",
    "HK1", "HK2", "PFKL", "PFKM", "ALDOA", "ALDOC",
    "PDHA1", "CS", "IDH1", "IDH2", "IDH3A", "SDHA", "SDHB", "FH", "MDH1", "MDH2",
    "ATP5F1A", "ATP5F1B", "ATP5F1C", "ATP5F1D", "ATP5F1E",
    "UQCRB", "UQCRC1", "UQCRC2", "COX4I1", "COX5A", "COX5B", "COX6A1", "COX6B1",
    "NDUFB3", "NDUFB4", "NDUFB6", "NDUFA1", "NDUFA4", "NDUFV1",
    # Mitochondrial (nuclear-encoded)
    "TOMM20", "TOMM40", "TIMM8A", "TIMM13", "SLC25A3", "SLC25A6", "VDAC1", "VDAC2",
    # Heat shock proteins / chaperones
    "HSPA1A", "HSPA1B", "HSPA5", "HSPA8", "HSPA9",
    "HSP90AA1", "HSP90AB1", "HSP90B1",
    "HSPD1", "HSPE1", "HSPH1",
    "DNAJB1", "DNAJC7", "STIP1", "BAG1", "BAG3",
    # Transcription factors (general)
    "JUN", "FOS", "FOSB", "JUNB", "JUND",
    "EGR1", "EGR2", "EGR3", "KLF4", "KLF6",
    "SP1", "SP3", "NFE2L2", "YBX1", "YBX3",
    # MAPK / ERK signaling
    "MAPK1", "MAPK3", "MAPK8", "MAPK14",
    "MAP2K1", "MAP2K2", "MAP3K1", "MAP3K7",
    "RAF1", "BRAF", "ARAF",
    "RPS6KA1", "RPS6KA3", "RPS6KB1",
    # PI3K / AKT / mTOR
    "AKT1", "AKT2", "AKT3",
    "PIK3CA", "PIK3CB", "PIK3CD", "PIK3R1", "PIK3R2",
    "PTEN", "MTOR", "RPTOR", "RICTOR",
    "TSC1", "TSC2", "RHEB", "S6K1",
    # JAK-STAT signaling (non-IFN components)
    "JAK1", "JAK2", "JAK3", "TYK2",
    "STAT2", "STAT3", "STAT5A", "STAT5B", "STAT6",
    "SOCS1", "SOCS2", "SOCS3",
    # NF-kB signaling
    "NFKB1", "NFKB2", "RELA", "RELB", "REL",
    "IKBKA", "IKBKB", "IKBKG", "CHUK",
    "TNFAIP3", "NFKBIA", "NFKBIB",
    # Apoptosis / BCL2
    "BCL2", "BCL2L1", "BCL2L11", "BAX", "BAK1", "BAD",
    "CASP3", "CASP7", "CASP8", "CASP9",
    "MCL1", "BID", "PMAIP1", "BBC3",
    # Cell cycle
    "CDK1", "CDK2", "CDK4", "CDK6",
    "CCNA2", "CCNB1", "CCND1", "CCND2", "CCND3", "CCNE1", "CCNE2",
    "CDKN1A", "CDKN1B", "CDKN2A", "CDKN2B",
    "RB1", "E2F1", "E2F2", "E2F3",
    # DNA damage / repair
    "TP53", "MDM2", "BRCA1", "BRCA2",
    "ATM", "ATR", "CHEK1", "CHEK2",
    "H2AX", "RAD51", "XRCC1", "LIG3",
    # Vesicle trafficking
    "RAB1A", "RAB5A", "RAB7A", "RAB11A", "RAB27A",
    "ARF1", "ARF3", "ARF6",
    "COPA", "COPB1", "COPB2",
    # Immune signaling (non-marker)
    "PTPN6", "PTPN11", "PTPN22",
    "INPP5D", "SHIP2",
    "LCK", "FYN", "ZAP70", "ITK",
    "BTK", "SYK", "BLNK",
    "PIK3CD", "PLCG1", "PLCG2",
    "RASGRP1", "RASGRP3",
    "LAT", "LAT2", "SLP76",
    # Cytokines and receptors (not markers)
    "IL1A", "IL1B", "IL6", "IL10", "IL12A", "IL12B", "IL15", "IL18",
    "TNF", "TNFRSF1A", "TNFRSF1B",
    "IFNAR1", "IFNAR2", "IFNGR1", "IFNGR2",
    "IL6R", "IL10RA", "IL10RB",
    "CCL2", "CCL3", "CCL4", "CCL5", "CCL7", "CCL19", "CCL20",
    "CXCL1", "CXCL2", "CXCL3", "CXCL8", "CXCL12",
    "CCR1", "CCR2", "CCR3", "CCR5", "CCR6", "CCR7",
    "CXCR1", "CXCR2", "CXCR3", "CXCR4", "CXCR5",
    # Adhesion / migration
    "ICAM1", "VCAM1", "MADCAM1",
    "ITGAL", "ITGB2", "ITGB1", "ITGAV",
    "SELP", "SELE",
    # Complement
    "C1QA", "C1QB", "C1QC", "C3", "C3AR1", "C5AR1",
    "CFB", "CFD", "CFH",
    # Antigen presentation
    "HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1", "HLA-DMA", "HLA-DMB",
    "CD74", "CIITA",
    # NK receptor / cytotoxicity
    "KLRC1", "KLRC2", "KLRC3", "KLRC4",
    "KLRK1", "NCR1", "NCR2", "NCR3",
    "PRF1", "GZMH", "GZMM",
    # Monocyte / macrophage specific
    "ITGB3", "CLEC4E", "CLEC7A", "MARCO",
    "MRC1", "CD163", "SIGLEC1",
    "TREM1", "TREM2",
    "LILRA1", "LILRA2", "LILRB1", "LILRB2", "LILRB4",
    # B cell specific
    "IGHG1", "IGHG2", "IGHA1", "IGKC", "IGLC2",
    "CD19", "CD38", "SLAMF7",
    "AICDA", "BCL6", "IRF4",
    # T cell activation / costimulation
    "CD28", "CD80", "CD86", "CTLA4",
    "CD40", "CD40LG",
    "ICOS", "ICOSLG",
    "TNFRSF9", "TNFRSF4",
    # Metabolism — lipids / cholesterol
    "LDLR", "HMGCR", "FASN", "ACACA", "ACLY",
    "SCD", "ELOVL1", "ELOVL2",
    "PPARA", "PPARG", "PPARD", "RXRA",
    # RNA metabolism
    "XPO1", "XPO5", "YTHDF1", "YTHDF2", "YTHDF3",
    "DDX5", "DDX17", "DDX21", "DHX9",
    "DCP1A", "DCP2", "XRN1",
    # Miscellaneous PBMC-expressed
    "ACTB", "GAPDH", "HMGB1", "HMGB2", "HMGN1", "HMGN2",
    "S100A4", "S100A6", "S100A10", "S100A11",
    "CALM1", "CALM2", "CALM3", "CALR",
    "STMN1", "CFL1", "PFN1",
    "HPRT1", "PPIA", "PPIB", "PPID",
    "TXN", "TXNRD1", "PRDX1", "PRDX2", "PRDX3",
    "SOD1", "SOD2", "CAT", "GPX1", "GPX4",
    "PCNA", "PARP1", "APEX1",
    "TOP1", "TOP2A", "TOP2B",
    "TERT", "DKC1", "NOP56",
    "NCL", "NPM1", "FBL",
    "LMNA", "LMNB1", "LMNB2",
    "H2BC3", "H3-3A", "H4C1",
]


# Build the full gene list
def build_gene_list():
    """Assemble a de-duplicated list of genes with real HGNC symbols for KEGG/GO compatibility."""
    all_marker_genes = []
    for genes in CELL_TYPE_MARKERS.values():
        all_marker_genes.extend(genes)

    final_genes = list(dict.fromkeys(all_marker_genes + IFN_DE_GENES + MT_GENES + BACKGROUND_GENES))
    return final_genes


# ---------------------------------------------------------------------------
# Count generation
# ---------------------------------------------------------------------------

def neg_binom_counts(mean, dispersion, n):
    """Sample from a negative binomial parameterised by mean and dispersion (1/r)."""
    r = 1.0 / dispersion
    p = r / (r + mean)
    return np.random.negative_binomial(r, p, size=n).astype(np.float32)


def generate_sample(
    genes: list,
    n_cells: int,
    cell_type_fracs: dict,
    ifn_fold_changes: dict,
    is_treatment: bool,
    per_sample_noise: float,
    rng: np.random.Generator,
) -> ad.AnnData:
    """
    Generate a single synthetic sample.

    Parameters
    ----------
    genes : list of str
        Gene names (length n_genes).
    n_cells : int
        Number of cells in the sample.
    cell_type_fracs : dict
        {cell_type: fraction} — must sum to 1.
    ifn_fold_changes : dict
        {gene: fold_change} — applied to treatment samples only.
    is_treatment : bool
        Whether this is a treatment (IFN-stimulated) sample.
    per_sample_noise : float
        Multiplicative noise scale between samples of the same condition.
    rng : np.random.Generator
        Seeded RNG for reproducibility.
    """
    n_genes = len(genes)
    gene_idx = {g: i for i, g in enumerate(genes)}

    # Assign cells to cell types
    cell_types = list(cell_type_fracs.keys())
    fracs = np.array([cell_type_fracs[ct] for ct in cell_types])
    fracs /= fracs.sum()
    n_per_type = rng.multinomial(n_cells, fracs)

    all_counts = np.zeros((n_cells, n_genes), dtype=np.float32)
    obs_cell_type = []
    cell_offset = 0

    for ct, n_ct in zip(cell_types, n_per_type):
        if n_ct == 0:
            continue
        markers = CELL_TYPE_MARKERS.get(ct, [])

        # Background expression for the 687-gene panel.
        # These genes are biologically relevant (ribosomal, metabolic, immune),
        # so most should be detected in each cell. A mean of 0.5-1.5 yields
        # ~50-70% genes detected per cell (~350-480 genes), well above the
        # pipeline default min_genes=200.
        base_mean = rng.uniform(0.5, 1.5, size=n_genes) * rng.uniform(
            1 - per_sample_noise, 1 + per_sample_noise
        )
        ct_counts = np.column_stack(
            [neg_binom_counts(m, 0.5, n_ct) for m in base_mean]
        )

        # Mitochondrial gene expression: higher mean and higher dispersion than
        # background, producing a realistic MT% range (~5-20%) across cells.
        for mt_gene in MT_GENES:
            if mt_gene in gene_idx:
                idx = gene_idx[mt_gene]
                mt_mean = rng.uniform(5.0, 9.0) * rng.uniform(
                    1 - per_sample_noise, 1 + per_sample_noise
                )
                ct_counts[:, idx] = neg_binom_counts(mt_mean, 1.5, n_ct)

        # Marker gene overexpression
        for marker in markers:
            if marker in gene_idx:
                idx = gene_idx[marker]
                marker_mean = rng.uniform(3.0, 12.0) * rng.uniform(
                    1 - per_sample_noise, 1 + per_sample_noise
                )
                ct_counts[:, idx] = neg_binom_counts(marker_mean, 0.3, n_ct)

        # IFN response overexpression in treatment
        if is_treatment:
            for gene, fc in ifn_fold_changes.items():
                if gene in gene_idx:
                    idx = gene_idx[gene]
                    base = ct_counts[:, idx].mean() + 0.5
                    # Apply to ~70% of cells (stochastic activation)
                    stimulated = rng.random(n_ct) < 0.70
                    ct_counts[stimulated, idx] = neg_binom_counts(
                        base * fc * rng.uniform(0.8, 1.2), 0.3, stimulated.sum()
                    )

        all_counts[cell_offset : cell_offset + n_ct] = ct_counts
        obs_cell_type.extend([ct] * n_ct)
        cell_offset += n_ct

    # Shuffle cell order
    perm = rng.permutation(n_cells)
    all_counts = all_counts[perm]
    obs_cell_type = [obs_cell_type[i] for i in perm]

    # Build AnnData
    X = sp.csr_matrix(all_counts)
    obs = pd.DataFrame(
        {"cell_type_sim": obs_cell_type},
        index=[f"cell_{i:05d}" for i in range(n_cells)],
    )
    var = pd.DataFrame(index=genes)
    var.index.name = "gene_symbols"

    adata = ad.AnnData(X=X, obs=obs, var=var)
    adata.var_names_make_unique()
    return adata


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Generate synthetic PBMC demo data for scAnnex (IFN-beta stimulation scenario)."
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default="data_demo",
        help="Output directory for generated H5AD files (default: data_demo/)",
    )
    parser.add_argument(
        "--n-cells",
        type=int,
        default=1000,
        help="Number of cells per sample (default: 1000)",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=42,
        help="Random seed for reproducibility (default: 42)",
    )
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    logger.info("Building gene list (~2000 human genes)...")
    genes = build_gene_list()
    logger.info(f"Gene panel: {len(genes)} genes ({len(IFN_DE_GENES)} IFN DE genes included)")

    # Cell type composition (PBMC-like proportions)
    cell_type_fracs = {
        "T_CD4": 0.25,
        "T_CD8": 0.18,
        "B_cell": 0.13,
        "NK_cell": 0.10,
        "Monocyte_CD14": 0.22,
        "Monocyte_CD16": 0.07,
        "other": 0.05,
    }

    # Fold changes for IFN DE genes (log2FC ~ 1.5-3.0 → FC ~ 3-8x)
    ifn_fcs = {}
    rng_fc = np.random.default_rng(args.seed)
    for gene in IFN_ALPHA_GENES:
        ifn_fcs[gene] = float(rng_fc.uniform(4.0, 10.0))
    for gene in IFN_GAMMA_GENES:
        if gene not in ifn_fcs:
            ifn_fcs[gene] = float(rng_fc.uniform(3.0, 8.0))

    # Sample definitions
    samples = [
        {"sample_id": "Ctrl_PBMC_C1", "batch": "batch1", "condition": "Control", "is_treatment": False, "seed_offset": 0},
        {"sample_id": "Ctrl_PBMC_C2", "batch": "batch2", "condition": "Control", "is_treatment": False, "seed_offset": 1},
        {"sample_id": "Ctrl_PBMC_C3", "batch": "batch3", "condition": "Control", "is_treatment": False, "seed_offset": 2},
        {"sample_id": "Treat_PBMC_T1", "batch": "batch1", "condition": "Treatment", "is_treatment": True, "seed_offset": 3},
        {"sample_id": "Treat_PBMC_T2", "batch": "batch2", "condition": "Treatment", "is_treatment": True, "seed_offset": 4},
        {"sample_id": "Treat_PBMC_T3", "batch": "batch3", "condition": "Treatment", "is_treatment": True, "seed_offset": 5},
    ]

    for meta in samples:
        sample_id = meta["sample_id"]
        rng = np.random.default_rng(args.seed + meta["seed_offset"] * 100)

        logger.info(f"Generating {sample_id} ({meta['condition']}, {meta['batch']}, {args.n_cells} cells)...")
        adata = generate_sample(
            genes=genes,
            n_cells=args.n_cells,
            cell_type_fracs=cell_type_fracs,
            ifn_fold_changes=ifn_fcs if meta["is_treatment"] else {},
            is_treatment=meta["is_treatment"],
            per_sample_noise=0.15,
            rng=rng,
        )

        # Attach sample-level metadata
        adata.obs["sample_id"] = sample_id
        adata.obs["batch"] = meta["batch"]
        adata.obs["condition"] = meta["condition"]

        # Rename cell barcodes to include sample ID (mimics CellRanger barcodes)
        adata.obs_names = [f"{sample_id}_{bc}" for bc in adata.obs_names]

        # Output filename convention from nf-core/scrnaseq mtx_conversions
        out_file = output_dir / f"{sample_id}_filtered_matrix.h5ad"
        adata.write_h5ad(out_file)
        logger.info(f"  Written: {out_file} ({adata.n_obs} cells x {adata.n_vars} genes)")

    logger.info("")
    logger.info("Demo data generation complete.")
    logger.info(f"Output directory: {output_dir.resolve()}")
    logger.info("Files generated:")
    for meta in samples:
        logger.info(f"  {meta['sample_id']}_filtered_matrix.h5ad  [{meta['condition']}, {meta['batch']}]")

    # Create logs directory so SLURM can write output files on first sbatch run
    logs_dir = output_dir / "logs"
    logs_dir.mkdir(exist_ok=True)

    logger.info("")
    logger.info("Next step: run the pipeline with")
    logger.info("  bash data_demo/run_command.sh")


if __name__ == "__main__":
    main()
