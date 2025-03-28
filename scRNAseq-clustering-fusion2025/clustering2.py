import scanpy as sc
import numpy as np
import pandas as pd
from sklearn.metrics import silhouette_score
from scipy.cluster.hierarchy import linkage, dendrogram
import itertools

# ... Load your AnnData object (adata) ...

adata= sc.read_h5ad("adata_filtered_combined_feb2025.h5ad")

#clustering starts here
n_iterations = 1000  # Increased iterations for robustness
length = len(adata.obs_names)
pairwise_cooccurrence = np.zeros((length, length))
indices = np.arange(length)
#print(pairwise_cooccurrence)

for i in range(n_iterations):
    # Subsample cells and genes (90%)
    cells_to_keep_indices = np.random.choice(indices, size=int(0.9 * length), replace=False)
    genes_to_keep = np.random.choice(adata.var_names, size=int(0.9 * len(adata.var_names)), replace=False)
    cells_to_keep = adata.obs_names[cells_to_keep_indices]
    adata_subsampled = adata[cells_to_keep, genes_to_keep].copy()

    # Leiden clustering (default parameters)
    sc.pp.neighbors(adata_subsampled)
    sc.tl.leiden(adata_subsampled)
    #cluster_counts =pd.group(adata_subsampled.obs["leiden"])
    total_cluster_num = adata_subsampled.obs["leiden"].nunique()
    #print(cluster_num)
    
    #expand leiden matrix back to normal dimensions and one-hot-encode
    leiden_matrix = np.zeros((length, total_cluster_num))
    for j in range(length):
        if j not in cells_to_keep_indices:
                continue
        cell = adata.obs_names[j]
        cluster = int(adata_subsampled.obs["leiden"][cell])
        leiden_matrix[j][cluster] = 1
    #print(leiden_matrix)
    
    # Track pairwise co-occurrence -- using one hot encoding!!!
    #print(test)
    #print(test.index)
    a = np.dot(leiden_matrix,leiden_matrix.T) #dot product counts number of times cells are paired together across ALL cells
    pairwise_cooccurrence = pairwise_cooccurrence + a
    #print(pairwise_cooccurrence)
    
    print("completed run: " + str(i), flush=True)
    #break

final_pairwise_df = pd.DataFrame(pairwise_cooccurrence)
final_pairwise_df.columns = adata.obs_names
final_pairwise_df.index = adata.obs_names

final_pairwise_df.to_csv("cell_pair_counts.csv")
