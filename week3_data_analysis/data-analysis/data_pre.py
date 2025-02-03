import h5py
import pandas as pd
from scipy.sparse import csc_matrix


def load_hdf5_as_dataframe(h5_file_path):
    """
    Loads an HDF5 file containing a sparse gene expression matrix and converts it to a dense pandas DataFrame.

    Parameters:
        h5_file_path (str): Path to the HDF5 file.

    Returns:
        pd.DataFrame: Gene expression matrix with rows as cell barcodes and columns as gene names.
    """
    with h5py.File(h5_file_path, "r") as f:
        # Extract matrix components
        matrix_group = f["matrix"]
        data = matrix_group["data"][:]
        indices = matrix_group["indices"][:]
        indptr = matrix_group["indptr"][:]
        shape = tuple(matrix_group["shape"][:])

        # Build sparse matrix
        sparse_matrix = csc_matrix((data, indices, indptr), shape=shape)

        # Extract barcodes (row identifiers) and gene names (column identifiers)
        cell_barcodes = matrix_group["barcodes"][:].astype(str)
        gene_names = matrix_group["features"]["name"][:].astype(str)

        # Convert sparse matrix to a dense DataFrame
        expression_df = pd.DataFrame(sparse_matrix.T.toarray(), index=cell_barcodes, columns=gene_names)

    return expression_df


def save_preprocessed_data(df, file_path):
    """
    Saves a gene expression DataFrame to a CSV file.

    Parameters:
        df (pd.DataFrame): Gene expression DataFrame to save.
        file_path (str): Path to save the CSV file.
    """
    df.to_csv(file_path)
    print(f"Data saved to {file_path}")


def prepare_data(v1_file_path, v2_file_path):
    """
    Loads and preprocesses data from two HDF5 files for Xenium v1 and Xenium v2.

    Parameters:
        v1_file_path (str): Path to the Xenium v1 HDF5 file.
        v2_file_path (str): Path to the Xenium v2 HDF5 file.

    Returns:
        tuple: Preprocessed DataFrames for Xenium v1 and Xenium v2.
    """
    print("Loading Xenium v1 data...")
    v1_data = load_hdf5_as_dataframe(v1_file_path)

    print("Loading Xenium v2 data...")
    v2_data = load_hdf5_as_dataframe(v2_file_path)

    return v1_data, v2_data


# Example Usage
if __name__ == "__main__":
    # File paths
    xenium_v1_file = "cell_feature_matrix_v1.h5"
    xenium_v2_file = "cell_feature_matrix_v2.h5"

    # Preprocess the data
    v1_matrix, v2_matrix = prepare_data(xenium_v1_file, xenium_v2_file)

    # Save preprocessed data
    save_preprocessed_data(v1_matrix, "xenium_v1_preprocessed.csv")
    save_preprocessed_data(v2_matrix, "xenium_v2_preprocessed.csv")