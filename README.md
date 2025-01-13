
code - https://github.com/qu-bit1/scRNAseq-analysis

data - https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE171524
## Structure of CSV 
The initial CSV file contains **raw gene expression counts** from a **single-cell RNA sequencing (scRNA-seq)** experiment.
- **Rows**: Represent **genes** (features).
- **Columns**: Represent individual **cells** (observations).

| Gene Name | Cell 1 | Cell 2 | Cell 3 | ... | Cell N |
| --------- | ------ | ------ | ------ | --- | ------ |
| GeneA     | 5      | 0      | 2      | ... | 10     |
| GeneB     | 0      | 3      | 1      | ... | 4      |
| GeneC     | 0      | 0      | 0      | ... | 0      |
| ...       | ...    | ...    | ...    | ... | ...    |
|           |        |        |        |     |        |
We transpose the matrix so that cells are represented as rows and genes as columns, which is the standard format for scRNA-seq analysis.

- **Gene Names**:  
    The first column contains the names of genes.
- **Cell Identifiers**:  
    Each column (after the first) represents a **single cell** captured during the sequencing process.
- **Gene Expression Counts**:  
    The values in the table are **raw counts** of gene expression, i.e., the number of times each gene’s transcript was observed in each cell. These counts are typically integers and represent the abundance of RNA molecules for each gene in each cell. Most of the values will be **zero**, as single cells only express a small subset of all possible genes at a given time, resulting in a **sparse matrix**.

Since the data contains **raw counts**:
1. Preprocessing steps like filtering lowly expressed genes and cells with abnormal characteristics are essential for removing noise.
2. Highly variable gene selection is critical because most genes will not contribute meaningful information to cell type differentiation.
3. Log normalization or other transformations are typically applied after the preprocessing steps in preparation for dimensionality reduction (e.g., PCA or UMAP) and clustering.

---
## AnnData obj
Highly efficient data structure designed for single-cell genomics data. We need it cause we need to store metadata too along with the count matrix.
- `adata.X` - count matrix 
- `adata.obs` - a df for **cell-level** metadata where each row corresponds to a cell in the dataset. Can include Sample IDs, QC metrics (e.g., number of genes detected, etc), Clustering results or cell type annotations.
- `adata.var` - for **gene-level** metadata
and many other like `uns` etc.

---
## Filtering Genes

1. **Min cell count** - Removes genes that are expressed in fewer than some 'n' cells. This step reduces noise by excluding rarely expressed genes that may not contribute meaningful information.
2. **Selecting highly variable genes** - Identifies the top 'n' most variable genes using the **Seurat v3** method. Highly variable genes are important for clustering and dimensionality reduction because they capture the most meaningful biological variation between cells.
---

## Setting Up and Training the SCVI Model**

```python
scvi.model.SCVI.setup_anndata(adata)
vae = scvi.model.SCVI(adata)
vae.train(accelerator='gpu', precision='16-mixed')
```

- **SCVI** (Single-Cell Variational Inference) is a probabilistic model that uses deep learning to represent scRNA-seq data in a lower-dimensional latent space.
- The model is trained using a **GPU** with **mixed precision** to speed up computation and take advantage of Tensor Cores trading off some accuracy for performance.
```python
import torch
torch.set_float32_matmul_precision('high')
```
- Can leave this out not necessary, just to speedup the A100 I was using, which didn't reached the capacity I wanted it to, probably some processing bottleneck in code.
- The latent representation produced by SCVI can later be used for clustering, visualization, or downstream analyses.

---

## Doublet Detection & Filtration with SOLO**

```python
solo = scvi.external.SOLO.from_scvi_model(vae)
solo.train(accelerator='gpu', precision='16-mixed')
df = solo.predict()
```

- **Doublets** are artificial cells created when two cells are captured together during sequencing. They can introduce false clusters and distort downstream analysis.
- **SOLO** is a doublet detection method based on SCVI. It trains a classifier on synthetic doublets and real singlets to distinguish between them.
- After training, predictions are made, and each cell is classified as either a "singlet" or a "doublet."
- Cells predicted as doublets with a confidence difference greater than 0.4 are selected for removal. 0.4 just out of pure interpretation from all the EDA.
---

## Annotating the Sample**

```python
adata.obs['Sample'] = csv_path.split('_')[2]
```

- Extracts a sample identifier from the file name and stores it in the observation metadata (`adata.obs`).

---

## Quality Control (QC) Filtering

### **Filtering by Minimum Genes per Cell**

```python
sc.pp.filter_cells(adata, min_genes=200)
```

- Removes cells that express fewer than 200 genes.
- This step helps exclude empty droplets or cells with insufficient RNA capture.

### **Annotating Mitochondrial and Ribosomal Genes**

```python
adata.var['mt'] = adata.var_names.str.startswith('mt-')  # Mitochondrial genes
adata.var['ribo'] = adata.var_names.isin(ribo_genes[0].values)  # Ribosomal genes
```

- Mitochondrial genes are identified by their prefix `mt-`.
- Ribosomal genes are identified using a list `ribo_genes`.

### **Calculating QC Metrics**

```python
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt', 'ribo'], percent_top=None, log1p=False, inplace=True)
```

- Computes QC metrics such as the percentage of mitochondrial and ribosomal gene expression for each cell.
- High mitochondrial content is often indicative of dying cells, and high ribosomal content can indicate technical artifacts.

### **Filtering by QC Thresholds**

```python
upper_lim = np.quantile(adata.obs.n_genes_by_counts.values, 0.98)
adata = adata[adata.obs.n_genes_by_counts < upper_lim]
adata = adata[adata.obs.pct_counts_mt < 20]
adata = adata[adata.obs.pct_counts_ribo < 2]
```

- Filters out cells with:
    - Extremely high gene counts (above the 98th percentile), which may represent multiplets.
    - Mitochondrial content greater than 20%, as high mitochondrial content suggests low-quality or dying cells.
    - Ribosomal content greater than 2%, as excessive ribosomal content can indicate technical noise.

Returning the Preprocessed Data. The cleaned, filtered AnnData object is returned, ready for downstream analyses such as clustering, dimensionality reduction, and visualization.

---