## Core Alignment Strategy

The alignment employs an affine transformation matrix that handles:
- Translation in x and y directions
- Rotation by theta angle
- Uniform scaling factor

The transformation is expressed mathematically as:
```
x_new = scale * (x * cos(θ) - y * sin(θ)) + tx
y_new = scale * (x * sin(θ) + y * cos(θ)) + ty
```

What makes this approach effective is its optimization method. The code uses **Powell's algorithm** to minimize the mean distance between transformed points and their nearest neighbors in the reference dataset. This avoids getting stuck in local minima that can plague gradient-based methods.

## Spatial Point Cloud Registration

The spatial alignment treats cell centroids as point clouds and optimizes their overlap. The key technical aspects:

1. The code builds a k-dimensional tree (cKDTree) for efficient spatial queries
2. It iteratively evaluates transformation parameters to minimize the objective function
3. The error function calculates mean distances between transformed points and their nearest matches

The initial parameters [0, 0, 0, 1.5] provide a starting point, with the scaling factor intentionally set to account for potential differences in coordinate scales.

## Gene Expression Integration

After spatial alignment, the code performs feature-based integration:

1. Identify common genes between datasets (intersection of feature spaces)
2. Preserve dataset origin information via batch annotations
3. Apply the optimized transformation to both cell centroids and individual transcript coordinates

This ensures that both cellular and subcellular data are properly aligned.
