import math
import time
import pandas as pd
import numpy as np
import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt
import paste as pst
import cupy as cp
from matplotlib import style


def load_slices(slice_names=["slice1", "slice2"], use_parquet=True):
    slices = []
    for slice_name in slice_names:
        if use_parquet:
            df = pd.read_parquet(f"{slice_name}.parquet")
            df = df.T
            slice_i = sc.AnnData(df)
        else:
            slice_i = sc.read_csv(f"{slice_name}.csv")

        slice_i_coor = np.genfromtxt(f"{slice_name}_coor.csv", delimiter=",")
        
        if slice_i_coor.shape[0] != slice_i.n_obs:
            raise ValueError(f"Spatial coordinates shape {slice_i_coor.shape} does not match the number of cells {slice_i.n_obs} in {slice_name}.")
        
        slice_i.obsm['spatial'] = slice_i_coor

        sc.pp.filter_genes(slice_i, min_counts=20)
        sc.pp.filter_cells(slice_i, min_counts=100)

        slices.append(slice_i)
    return slices


def align_tissues(slice1, slice2, device='cpu'):
    # Convert to GPU arrays if device is GPU
    if device == 'gpu':
        slice1_X = cp.array(slice1.X)
        slice2_X = cp.array(slice2.X)
    else:
        slice1_X = np.array(slice1.X)
        slice2_X = np.array(slice2.X)

    result = pst.pairwise_align(slice1_X, slice2_X)  # Perform pairwise alignment

    result_centered = result - cp.mean(result, axis=0)

    return result_centered

# Function to plot and save alignment images
def plot_and_save_alignment(slice1, slice2, alignment_result, slice_colors, output_dir="output/"):
    fig, axs = plt.subplots(2, 2, figsize=(7,7))

    # Plot aligned slice 1
    pst.plot_slice(slice1, slice_colors[0], ax=axs[0,0])
    axs[0,0].set_title("Slice 1 Aligned")

    # Plot aligned slice 2
    pst.plot_slice(slice2, slice_colors[1], ax=axs[0,1])
    axs[0,1].set_title("Slice 2 Aligned")

    # Plot the alignment result
    sns.heatmap(cp.asnumpy(alignment_result), ax=axs[1,0], cmap="viridis")
    axs[1,0].set_title("Alignment Result")

    # Save the figure
    plt.savefig(f"{output_dir}/slice_alignment.png", dpi=300, bbox_inches='tight')
    plt.close(fig)


# Main execution logic
def main(slice_names=["slice1", "slice2"], use_parquet=True, device='gpu', output_dir="output/"):
    print("Loading slices...")
    slices = load_slices(slice_names=slice_names, use_parquet=use_parquet)
    slice1, slice2 = slices

    print("Aligning slices...")
    aligned_result = align_tissues(slice1, slice2, device=device)

    print("Saving plots...")
    plot_and_save_alignment(slice1, slice2, aligned_result, slice_colors=['#e41a1c','#377eb8'], output_dir=output_dir)
    print("Alignment completed and saved!")

if __name__ == "__main__":
    main(slice_names=["slice1", "slice2"], use_parquet=True, device='gpu', output_dir="output/")
