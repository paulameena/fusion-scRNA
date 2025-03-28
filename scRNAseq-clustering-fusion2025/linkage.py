import scanpy as sc
import numpy as np
import pandas as pd
from sklearn.metrics import silhouette_score
from scipy.cluster.hierarchy import linkage, dendrogram
import itertools

#load clustering output
cell_pair_counts = pd.read_csv("/mnt/pan/SOM_CCCC_JGS25/shultesp/scRNAseq-clustering-fusion2025/cell_pair_counts.csv", index_col=0)

#Part2. Frequency Calculation and Hierarchical Clustering:

# Convert to a distance matrix (1-normalized co-occurrence)
cooccurrence_matrix_normalized = cell_pair_counts.div(1000) #divide by number of iterations
distances = 1 - cooccurrence_matrix_normalized
#distances

# Hierarchical clustering using ward linkage (can experiment with other methods)
linkage_matrix = linkage(distances, method='ward')

# Convert the array to a Pandas DataFrame
df = pd.DataFrame(linkage_matrix)

# Save the DataFrame to a CSV file
csv_file_path = 'linkage.csv'
df.to_csv(csv_file_path, index=False)
