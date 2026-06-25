#!/public/home/gzcui_gdl/miniconda3/envs/r-reticulate/bin/python
###############################################################################
# Feature Selection Pipeline for High-Dimensional Data
# 
# This script performs multiple feature selection algorithms on input data
# to identify the most informative features (e.g., genes). It integrates
# 5 different feature selection methods from the custom feature_selection module.
#
# Usage: python feature.selection.byhand.py <input_matrix.csv> <num_features> <num_classes>
#
# Args:
#   input_matrix.csv: Path to CSV file with features as rows and samples as columns
#   num_features: Number of top features to select (nfs)
#   num_classes: Number of classes/clusters for supervised methods (ncls)
#
# Output:
#   Five CSV files containing selected feature indices/rankings from each method:
#   - fs_lap_score.<nfs>.csv   : Laplacian Score-based selection
#   - fs_low_variance.<nfs>.csv: Low Variance-based selection
#   - fs_SPEC.<nfs>.csv        : SPEC-based selection
#   - fs_MCFS.<nfs>_<ncls>.csv : Multi-Cluster Feature Selection
#   - fs_NDFS.<nfs>_<ncls>.csv : Neighborhood Discriminant Feature Selection
#

###############################################################################

import pandas as pd
import numpy as np
import feature_selection as fs  # Custom feature selection module
import time
import os
import sys

# Parse command-line arguments
tmat = sys.argv[1]       # Input matrix file path (CSV format)
nfs = int(sys.argv[2])   # Number of top features to select (nfs)
ncls = int(sys.argv[3])  # Number of classes/clusters for supervised methods (ncls)
print("Input matrix:", tmat)

# Load input data from CSV file
# Expected format: rows = features (genes), columns = samples
df = pd.read_csv(tmat, index_col=0)

# Print data preview for verification
print("Data preview:")
print(df.head())

###############################################################################
# Feature Selection Algorithms
# 
# Execute 5 different feature selection methods sequentially.
# Each method returns the indices/rankings of selected features.
###############################################################################

# 1. Laplacian Score (fs_lap_score)
# Unsupervised method that measures the locality-preserving power of features
# by computing the Laplacian score based on graph Laplacian
data1 = fs.fs_lap_score(df.values, nfs)
print(time.strftime("%Y-%m-%d %H:%M:%S") + ": fs_lap_score completed")

# 2. Low Variance Filter (fs_low_variance)
# Unsupervised method that removes features with low variance
# Features with variance below a threshold are considered uninformative
data2 = fs.fs_low_variance(df.values, nfs)
print(time.strftime("%Y-%m-%d %H:%M:%S") + ": fs_low_variance completed")

# 3. SPEC (fs_SPEC)
# Spectral feature selection method that uses the eigenvalues and eigenvectors
# of the similarity matrix to select features that best preserve the global structure
data3 = fs.fs_SPEC(df.values, nfs)
print(time.strftime("%Y-%m-%d %H:%M:%S") + ": fs_SPEC completed")

# 4. Multi-Cluster Feature Selection (fs_MCFS)
# Semi-supervised method that considers the cluster structure of data
# Selects features that are discriminative across multiple clusters
# Requires ncls (number of classes) parameter
data4 = fs.fs_MCFS(df.values, nfs, ncls)
print(time.strftime("%Y-%m-%d %H:%M:%S") + ": fs_MCFS completed")

# 5. Neighborhood Discriminant Feature Selection (fs_NDFS)
# Semi-supervised method that maximizes the margin between different classes
# while preserving the local structure of data
# Requires ncls (number of classes) parameter
data5 = fs.fs_NDFS(df.values, nfs, ncls)
print(time.strftime("%Y-%m-%d %H:%M:%S") + ": fs_NDFS completed")

###############################################################################
# Save Results
# 
# Save selected feature indices/rankings to CSV files for downstream analysis
###############################################################################

# Convert parameters to strings for filename construction
Nfs = str(nfs)
Ncls = str(ncls)

# Save results from each method
np.savetxt('fs_lap_score.' + Nfs + '.csv', data1)          # Laplacian Score results
np.savetxt('fs_low_variance.' + Nfs + '.csv', data2)       # Low Variance results
np.savetxt('fs_SPEC.' + Nfs + '.csv', data3)               # SPEC results
np.savetxt('fs_MCFS.' + Nfs + '_' + Ncls + '.csv', data4)  # MCFS results (includes ncls)
np.savetxt('fs_NDFS.' + Nfs + '_' + Ncls + '.csv', data5)  # NDFS results (includes ncls)



