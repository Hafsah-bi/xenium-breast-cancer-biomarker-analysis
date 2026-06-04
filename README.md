# xenium-breast-cancer-biomarker-analysis
Reproduces Xenium spatial transcriptomics analyses for breast cancer biomarker quantification, including housekeeping gene stability, LDHA normalization, myoepithelial markers, and periductal MMP11 analysis.

Reproduction of Figures 1A, 1B, 1C, 1E, 1F, 2C, 2F, 2G from:

> Janesick A, et al. *Biomarker Quantification in Breast Cancer using Xenium In Situ.* bioRxiv (2025). https://doi.org/10.64898/2025.12.08.692193

**Dataset:** [10x Genomics — Xenium FFPE Human Breast Biomarkers](https://www.10xgenomics.com/datasets/xenium-ffpe-human-breast-biomarkers)

---

## Table of Contents

1. [Overview](#1-overview)
2. [Project Structure](#2-project-structure)
3. [Input Data](#3-input-data)
4. [Setup & Installation](#4-setup--installation)
5. [Code Walkthrough](#5-code-walkthrough)
6. [Outputs](#6-outputs)
7. [Key Tables](#7-key-tables)
8. [Dependencies](#8-dependencies)

---

## 1. Overview

### 1.1 Purpose

Two goals:

1. **Identify stable housekeeping (HK) genes** at single-cell resolution — genes whose expression is consistent across cell types and tissue sections, suitable as normalization references.
2. **Identify biomarkers that track DCIS progression** — genes that vary predictably as Ductal Carcinoma In Situ advances from low- to high-grade.

The platform — **Xenium In Situ** — measures individual RNA transcripts at their exact spatial location within FFPE breast tissue, at single-cell resolution. This enables cell-type-specific analysis and spatial mapping of gene expression relative to duct boundaries.

---

### 1.2 Background

**Why DCIS matters:**  
DCIS is an early-stage breast cancer confined within milk ducts. Roughly 20–50% of untreated cases may progress to invasive cancer, but no reliable molecular predictor currently exists. Better biomarkers could reduce over-treatment.

**The normalization problem:**  
Spatial transcriptomics measurements are affected by tissue preparation quality, RNA degradation from FFPE fixation, and capture efficiency variation between samples. Normalization against stable reference genes corrects for these technical factors — but only if those reference genes are actually stable.

**Why bulk-derived HK genes fail at single-cell resolution:**  
Genes like GAPDH and TFRC were selected as Oncotype Dx references using bulk RNA data, where signals from many cells are averaged. At single-cell resolution, GAPDH expression rises with tumor aggressiveness (Warburg effect), making it unsuitable as a reference — normalizing against it can reverse the apparent trend of genuine biomarkers like LDHA.

---

### 1.3 Methodology

**Study design:** 12 FFPE human breast tissue sections (8 DCIS, 3 invasive ductal carcinoma, 1 normal), profiled with a custom 280-gene Xenium panel including HK candidates, cell-type markers, clinical biomarkers (ESR1, PGR, ERBB2, MKI67), and basement membrane genes.

**Pipeline summary:**

| Step | Description |
|------|-------------|
| Cell segmentation | Xenium on-instrument pipeline v4.0.0 → cell boundaries + per-cell transcript counts |
| QC filtering | Retain cells with >10 total transcripts and >5 unique genes |
| Normalization | SCTransform → PCA → neighbor graph → UMAP → Leiden clustering (resolution 0.3) |
| HK stability | CV (std/mean) computed intra-section (per cell type) and inter-section (across 12 sections) |
| Biomarker validation | LDHA normalization compared across HK gene sets |
| Myoepithelial analysis | LAMC2, KRT14, COL17A1 in normal vs. DCIS-adjacent myoepithelial cells |
| Spatial MMP11 | Cell-agnostic transcript density in 30 µm periductal rings around manually annotated ROIs |

**HK gene selection outcome:**  
Four genes emerged as best references for tumor cells: **EEF1G**, **EEF2**, **MALAT1**, **RPLP0** — all showing low intra- and inter-section CV across cell types. GAPDH, GUSB, and TFRC were found unsuitable.

**MMP11 spatial analysis:**  
Many stromal transcripts go unassigned to cells due to poor segmentation of thin, irregular fibroblasts in the ECM. To capture these, the authors used a cell-agnostic approach: manually delineated 16 duct ROIs in Xenium Explorer (8 MMP11-high, 8 MMP11-low), built 30 µm periductal rings around each, and counted all MMP11 transcripts inside regardless of cell assignment.

---

### 1.4 Key Findings

| Finding | Gene(s) | Notes |
|---------|---------|-------|
| Stable HK genes (single-cell) | EEF1G, EEF2, MALAT1, RPLP0 | Low CV intra- and inter-section |
| Unstable Oncotype Dx HK genes | GAPDH, GUSB, TFRC | Warburg-affected or high variability |
| Grade-increasing biomarkers | LDHA, SDC1 | Rise with DCIS progression |
| Grade-decreasing biomarkers | SFRP1, PIGR | Fall with DCIS grade |
| Invasion marker (myoepithelial) | LAMC2 | Upregulated in DCIS-adjacent myoepithelium |
| Peritumoral invasion signal | MMP11 | High periductal density correlates with MKI67+, CCNB1+, and aggressive signatures |

---

## 2. Project Structure

```
xenium_breast_cancer/
│
├── Xenium_Breast_Cancer_Biomarker_Analysis.ipynb
│
├── xenium_data/                        # Input (download separately)
│   ├── cell_feature_matrix.zarr.zip    # 25.1 MB
│   ├── cells.zarr.zip                  # 206.2 MB
│   ├── transcripts.zarr.zip            # 1.18 GB
│   ├── gene_panel.json                 # 153 KB
│   └── table.xlsx                      # Table S5 — ROI polygon coordinates
│
└── xenium_outputs/                     # Generated by notebook
    ├── qc_distributions.png
    ├── fig1A_hk_cv_violin.png
    ├── fig1B_deg_cv_violin.png
    ├── fig1C_cv_scatter.png
    ├── fig1E_LDHA_raw.png
    ├── fig1F_LDHA_normalized.png
    ├── fig2C_myoepithelial_markers.png
    ├── fig2F_MMP11_raw_density.png
    ├── fig2G_MMP11_normalized_density.png
    ├── xenium_annotated.h5ad
    ├── hk_gene_cv_table.csv
    └── tableS3_stable_genes.csv
```

---

## 3. Input Data

Download from the [10x Genomics dataset page](https://www.10xgenomics.com/datasets/xenium-ffpe-human-breast-biomarkers) — specifically the **S1-Top** section:

```
https://cf.10xgenomics.com/samples/xenium/4.0.0/Human_Breast_Biomarkers_S1_Top/
Human_Breast_Biomarkers_S1_Top_xe_outs.zip
```

| File | Size | Format | Contents |
|------|------|--------|----------|
| `cell_feature_matrix.zarr.zip` | 25.1 MB | Zarr | Sparse CSC matrix: 209,467 cells × 542 features |
| `cells.zarr.zip` | 206.2 MB | Zarr | Per-cell metadata: centroid (x, y), area, nucleus data |
| `transcripts.zarr.zip` | 1.18 GB | Zarr | Individual transcript records: coordinates, gene, cell assignment |
| `gene_panel.json` | 153 KB | JSON | Panel definition: 280+ targets with probe metadata |
| `table.xlsx` | — | Excel | Table S5: 16 ROI polygon coordinates for MMP11 analysis |

---

## 4. Setup & Installation

### 4.1 Google Drive

Upload all five files to `xenium_data/` in your Google Drive:

```
My Drive/
└── xenium_data/
    ├── cell_feature_matrix.zarr.zip
    ├── cells.zarr.zip
    ├── transcripts.zarr.zip
    ├── gene_panel.json
    └── table.xlsx
```

### 4.2 Mount in Colab

```python
from google.colab import drive
drive.mount('/content/drive')
```

### 4.3 Install Dependencies

```bash
!pip install -q zarr==2.18.3 scanpy anndata shapely geopandas scipy \
             matplotlib seaborn pandas numpy tqdm openpyxl
!pip install igraph leidenalg
```

### 4.4 Configure Paths

```python
DATA_DIR = '/content/drive/MyDrive/xenium_data'   # update if needed
LOCAL    = '/content/xenium_local'                 # temp extraction
OUT_DIR  = '/content/xenium_outputs'               # figure output
```

---

## 5. Code Walkthrough

### Section 0 — Install Dependencies
Installs all required packages. Run first in every new Colab session.

---

### Section 1 — Imports & Config

Defines the gene sets used throughout:

| Variable | Genes | Purpose |
|----------|-------|---------|
| `HK_PAPER` | EEF1G, EEF2, MALAT1, RPLP0 | Paper's recommended HK genes |
| `HK_ONCOTYPE` | ACTB, GAPDH, GUSB, RPLP0, TFRC | Oncotype Dx HK genes (comparison) |
| `HK_STABLE` | 15 genes | Full candidate set for Fig 1A |
| `DEG_GENES` | KIT, TFPI2, PIP, CHI3L1, etc. | High-CV DEGs for Fig 1B |
| `BIOMARKERS` | LDHA, SDC1, SFRP1, PIGR | Tumor grade biomarkers for Figs 1E/1F |
| `PIXEL_PER_UM` | 4.70588 | Xenium pixel size: 0.2125 µm/px |

---

### Section 2 — Mount Drive & Extract Zarr Archives
Mounts Drive and extracts `.zarr.zip` archives to `/content/xenium_local/`. Archives are extracted only once; skips if already present.

---

### Section 3 — Load Cell × Gene Matrix → AnnData

Builds the AnnData from the CSC zarr matrix:
1. Gene names from `gene_panel.json`
2. `data`, `indices`, `indptr` loaded from zarr → `scipy.sparse.csc_matrix` (genes × cells), then transposed
3. Cell IDs from two-column uint32 `cell_id` array
4. Cell metadata (centroids, areas) attached as `obs`
5. Spatial coordinates stored in `obsm['spatial']`

**Output:**
```
AnnData: 209,467 cells × 542 genes
obsm: 'spatial'  ← (x, y) centroid coordinates
obs:  cell_centroid_x/y, cell_area, nucleus_area
```

---

### Section 4 — Load Cell Metadata
Loads additional metadata from `cells.zarr`. Falls back to positional alignment if cell IDs between stores don't match.

---

### Section 5 — QC Filtering

Thresholds from paper:
- Keep cells with **> 10 total transcripts**
- Keep cells with **> 5 unique genes**

```
Before: 209,467 cells
After:  ~175,000–195,000 cells
```

Rebuilds AnnData with full 542-feature count, padding unnamed features as `blank_N`. Saves transcript and gene count distribution plots.

---

### Section 6 — Preprocessing Pipeline

| Step | Function | Settings |
|------|----------|---------|
| Preserve raw counts | `adata.layers['raw']` | Used for CV calculations |
| Normalize | `sc.pp.normalize_total` | Target sum = 10,000 |
| Log-transform | `sc.pp.log1p` | Stabilizes variance |
| HVG | `sc.pp.highly_variable_genes` | Top 2,000; flavor = `seurat` |
| PCA | `sc.pp.pca` | 30 components |
| Neighbors | `sc.pp.neighbors` | k-NN graph |
| UMAP | `sc.tl.umap` | 2D embedding |
| Clustering | `sc.tl.leiden` | Resolution = 0.3 |

---

### Section 7 — Cell Type Annotation

Scores each cell against 7 marker sets via `sc.tl.score_genes`; assigns the highest-scoring label:

| Cell Type | Marker Genes |
|-----------|-------------|
| Tumor | EPCAM, KRT8, KRT18, MKI67, KRT19 |
| Myoepithelial | KRT14, KRT5, ACTA2, MYLK, TP63 |
| T cell | CD3D, CD3E, CD8A, CD4, IL7R |
| B cell | CD79A, MS4A1, CD19 |
| Macrophage | CD68, CSF1R, LYZ |
| Fibroblast | COL1A1, DCN, FAP, VIM |
| Endothelial | PECAM1, VWF, CDH5 |

**Cell type counts (representative):**

| Cell Type | Count |
|-----------|------:|
| Tumor | 63,183 |
| Myoepithelial | 35,142 |
| Macrophage | 34,125 |
| Endothelial | 32,811 |
| T cell | 19,552 |
| B cell | 19,200 |

**Section 7b — Tumor Subtypes:**  
The paper uses 4 matched sections (S2T, S3T, S1B, S2B). This notebook sub-clusters tumor cells at resolution 0.5 and assigns subtype labels cyclically. This is an approximation — true subtype assignment requires matched multi-section data.

---

### Section 8 — CV Computation Helpers

**`get_raw_counts(adata, layer='raw')`**  
Returns the raw count matrix as a dense NumPy array; handles both sparse and dense formats.

**`compute_cv_per_celltype(adata, genes, ...)`**  
For each gene × cell type pair (minimum 5 cells), computes `CV = std / mean` from raw counts. Returns a `(genes × cell_types)` DataFrame.

---

### Section 9 — Figure 1A: HK Gene CV Violin Plot

CV distribution for 15 candidate HK genes across cell types, shown as violins with individual data points.

Color coding:
- Blue: Stable HK genes (paper recommendation)
- Green: Stable in tumor, not classified as universal HK
- Orange: Oncotype Dx HK genes

EEF1G, EEF2, MALAT1 should show tight low-CV violins; TFRC and GUSB should be wider and higher.

**Output:** `fig1A_hk_cv_violin.png`

---

### Section 10 — Figure 1B: DEG CV Violin Plot

CV of 9 differentially expressed genes across tumor subtypes. High CV here is expected — these genes are *supposed* to vary. All violins in magenta.

**Output:** `fig1B_deg_cv_violin.png`

---

### Section 11 — Figure 1C: CV vs Mean Scatter (Tumor Cells)

Log-log scatter of CV vs. mean expression for all 542 genes in tumor cells. Dashed crosshairs mark the stable, high-expression zone. Paper HK genes (EEF1G, EEF2, MALAT1, RPLP0) cluster in the low-CV, high-expression quadrant.

**Output:** `fig1C_cv_scatter.png`

---

### Section 12 — Figure 1E: LDHA Raw Expression

Raw LDHA transcript counts per tumor cell across subtypes (S2T → S3T → S1B → S2B). Jittered scatter with mean ± SEM and Mann-Whitney U significance brackets. Expected: progressive increase with grade (Warburg effect).

**Output:** `fig1E_LDHA_raw.png`

---

### Section 13 — Figure 1F: LDHA Normalized to HK Genes

Four-panel figure: LDHA normalized to EEF1G, RPLP0, GAPDH, and GUSB respectively. The critical comparison:
- EEF1G and RPLP0 normalization → preserves the grade-related LDHA trend
- GAPDH normalization → may flatten or reverse the trend (GAPDH itself rises with Warburg)
- GUSB normalization → may also distort the trend

**Output:** `fig1F_LDHA_normalized.png`

---

### Section 14 — Figure 2C: Myoepithelial Marker Expression

Bar chart of KRT14, COL17A1, and LAMC2 mean expression in myoepithelial cells adjacent to normal vs. DCIS ducts. Sub-clusters myoepithelial cells and assigns the first half as "Normal", second half as "Tumor-associated" (approximate; single-section limitation). LAMC2 label shown in red — the key upregulated invasion-associated marker.

**Output:** `fig2C_myoepithelial_markers.png`

---

### Section 15 — Figures 2F & 2G: MMP11 Periductal Transcript Density

Most complex section. Requires `table.xlsx` (Table S5) in `DATA_DIR`.

**Pipeline:**

```
Table S5 ROI polygons (16 total: 8 MMP11+, 8 MMP11−)
         ↓
Build Shapely Polygon per ROI
         ↓
Expand outward 30 µm (= 30 × 4.706 = 141.2 px)
Subtract inner duct → periductal ring
         ↓
Load transcripts.zarr (MMP11 + HK genes)
         ↓
Per ring:
  → Spatial join: find all transcripts within ring
  → Exclude transcripts assigned to Tumor or Myoepithelial cells
  → Count MMP11 transcripts
  → Count HK transcripts (EEF1G, EEF2, MALAT1, RPLP0)
         ↓
Fig 2F: raw_density    = MMP11_count / ring_area_µm²
Fig 2G: norm_density   = MMP11_count / geomean(HK_counts) / ring_area_µm²
```

**Pixel conversion:**  
Xenium pixel size = 0.2125 µm/px → 4.706 px/µm → 30 µm ring = 141.2 px outward expansion.

**Statistical test:** Mann-Whitney U — MMP11+ vs MMP11− rings.

**Outputs:** `fig2F_MMP11_raw_density.png`, `fig2G_MMP11_normalized_density.png`

---

### Section 16 — Save Outputs

| File | Contents |
|------|----------|
| `xenium_annotated.h5ad` | AnnData with cell type labels, UMAP, Leiden clusters |
| `hk_gene_cv_table.csv` | CV per gene per cell type (15 HK candidates) |
| `tableS3_stable_genes.csv` | All genes with min CV < 1.5, sorted |

---

## 6. Outputs

| File | Figure | Description |
|------|--------|-------------|
| `qc_distributions.png` | — | Transcript and gene count histograms post-filter |
| `fig1A_hk_cv_violin.png` | Fig 1A | CV of 15 HK candidate genes across cell types |
| `fig1B_deg_cv_violin.png` | Fig 1B | CV of 9 DEGs across tumor subtypes |
| `fig1C_cv_scatter.png` | Fig 1C | CV vs. mean expression (all genes, tumor cells) |
| `fig1E_LDHA_raw.png` | Fig 1E | Raw LDHA Tx/cell across 4 tumor subtypes |
| `fig1F_LDHA_normalized.png` | Fig 1F | LDHA normalized to EEF1G, RPLP0, GAPDH, GUSB |
| `fig2C_myoepithelial_markers.png` | Fig 2C | KRT14, COL17A1, LAMC2 in normal vs. DCIS-adjacent myoepithelium |
| `fig2F_MMP11_raw_density.png` | Fig 2F | MMP11 Tx/µm² in 30 µm periductal ring |
| `fig2G_MMP11_normalized_density.png` | Fig 2G | HK-normalized MMP11 density in periductal ring |
| `xenium_annotated.h5ad` | — | Annotated AnnData (reloadable) |
| `hk_gene_cv_table.csv` | — | Gene × cell-type CV table |
| `tableS3_stable_genes.csv` | — | Stable gene list (CV < 1.5) |

---

## 7. Key Tables

### Cell Count Summary

| Stage | Count |
|-------|------:|
| Raw cells in matrix | 209,467 |
| After QC | ~175,000–195,000 |
| Leiden clusters (res = 0.3) | ~8–15 |

### HK Gene Selection Summary

| Gene | Category | Recommended |
|------|----------|-------------|
| EEF1G | Translation elongation factor | Yes — very low CV |
| EEF2 | Translation elongation factor | Yes — very low CV |
| MALAT1 | lncRNA | Yes — very low CV |
| RPLP0 | Ribosomal protein | Yes — low CV |
| ACTB | Cytoskeletal | Borderline |
| GAPDH | Glycolytic enzyme | No — Warburg-affected |
| GUSB | Lysosomal enzyme | No — variable |
| TFRC | Transferrin receptor | No — high variability |

### Pixel-to-Micron Conversion

| Quantity | Value |
|----------|-------|
| Xenium pixel size | 0.2125 µm/px |
| Pixels per micron | 4.706 px/µm |
| 30 µm ring expansion | 141.2 px |
| Typical ring area | 50,000–200,000 µm² |

---

## 8. Dependencies

| Package | Version | Purpose |
|---------|---------|---------|
| `zarr` | 2.18.3 (pinned) | Read Xenium zarr archives |
| `scanpy` | ≥1.9 | Single-cell analysis |
| `anndata` | ≥0.9 | AnnData container |
| `numpy` | ≥1.24 | Array operations |
| `pandas` | ≥1.5 | DataFrames |
| `scipy` | ≥1.9 | Sparse matrices, statistics |
| `matplotlib` | ≥3.6 | Plotting |
| `seaborn` | ≥0.12 | Statistical visualization |
| `shapely` | ≥2.0 | Polygon geometry (periductal rings) |
| `geopandas` | ≥0.13 | Spatial join of transcripts to rings |
| `tqdm` | any | Progress bars |
| `igraph` | any | Graph backend for Leiden |
| `leidenalg` | any | Leiden clustering |
| `openpyxl` | any | Read Table S5 Excel file |
