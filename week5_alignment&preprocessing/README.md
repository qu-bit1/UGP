## Mapping/Aligning cells from both the datasets by Transcript Density-Based Mapping (Spatial Proximity + Expression Similarity)

Since two datasets may have different cell segmentation results, direct cell-to-cell matching is difficult

Map cells between two transcriptomic datasets by leveraging both:
1.	Spatial proximity (transcripts are physically located in space).
2.	Expression similarity (genes expressed in similar patterns should belong to the same cell).
Transcript Density-Based Mapping (Spatial Proximity + Expression Similarity)


**Step 1: Compute Transcript “Centers”**

Each transcript (dfT1, dfT2) has an associated (x, y) coordinate and belongs to a cell (cell_id).
However, since the cell segmentation may differ between the two datasets, we compute a “**center of mass**” for each cell.

This is done as:

$$x_{\text{center}} = \frac{\sum x_{\text{transcripts}}}{N_{\text{transcripts}}}, \quad y_{\text{center}} = \frac{\sum y_{\text{transcripts}}}{N_{\text{transcripts}}}$$

This step converts sparse transcript data into a more structured spatial representation of cells.

**Step 2: Find Closest Cells Between Datasets**

For each cell in Dataset 1 (dfC1), find the closest cell in Dataset 2 (dfC2) based on Euclidean distance.


## Further exploration
For alignment we can explore
**Transcript Profile Similarity (Expression Overlap + Correlation)**
- Each cell has a set of transcripts (genes), which define its molecular identity.
- We can create a gene expression profile for each cell and compute similarity using cosine or pearson correlation.

or a hybrid approach where we combine spatially and then refine by profile similarity