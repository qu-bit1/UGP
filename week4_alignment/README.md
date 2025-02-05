## Xenium File Alignment  

We employed two approaches for aligning Xenium files:  

### 1. **STalign**  
**Reference**: [STalign: Spatial Transcriptomics Alignment](https://www.nature.com/articles/s41467-023-43915-7)  
- A computational method for spatial alignment of transcriptomics data.  
- Based on Large Deformation Diffeomorphic Metric Mapping (LDDMM) and image varifolds.

### 2. **PASTE**  
**Reference**: [PASTE: Probabilistic Alignment of Spatial Transcriptomics Experiments](https://www.nature.com/articles/s41592-022-01459-6)   
- Leverages both transcriptional similarity and spatial distances between spots.  
- Computes pairwise slice alignment across adjacent ST slices.  
- Performs center slice integration by finding a low-rank expression matrix and mapping spots to other slices.  

