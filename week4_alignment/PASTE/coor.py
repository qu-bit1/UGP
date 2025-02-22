import pandas as pd
import numpy as np

cells_df = pd.read_parquet("cells.parquet")

print(cells_df.columns)

if "x_centroid" in cells_df.columns and "y_centroid" in cells_df.columns:
    coordinates = cells_df[["x_centroid", "y_centroid"]].values
    np.savetxt("slice2_coor.csv", coordinates, delimiter=",")
else:
    print("Columns 'x_centroid' and 'y_centroid' not found in cells.parquet!")
