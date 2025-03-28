#! /home/pvs13/.usr/local/python/3.10.4

import scanpy as sc
import numpy as np
import pandas as pd
from sklearn.metrics import silhouette_score
from scipy.cluster.hierarchy import linkage, dendrogram
import itertools


# ... Load your AnnData object (adata) ...

adata= sc.read_h5ad("adata_filtered_combined_feb2025.h5ad")

#run bootstrapping + repeated leiden clustering
n_iterations = 1000  # Increased iterations for robustness
pairwise_cooccurrence = np.zeros((len(adata.obs_names),len(adata.obs_names)))
original_list_indices=adata.obs_names.to_list()
#print(pairwise_cooccurrence)
#print(len(adata.obs_names))

for i in range(n_iterations):
    # Subsample cells and genes (90%)
    cells_to_keep = np.random.choice(adata.obs_names, size=int(0.9 * len(adata.obs_names)), replace=False)
    genes_to_keep = np.random.choice(adata.var_names, size=int(0.9 * len(adata.var_names)), replace=False)
    adata_subsampled = adata[cells_to_keep, genes_to_keep].copy()

    # Leiden clustering (default parameters)
    sc.pp.neighbors(adata_subsampled)
    sc.tl.leiden(adata_subsampled)

    # Track pairwise co-occurrence -- using one hot encoding!!!
    test = pd.get_dummies(adata_subsampled.obs['leiden'])
    #print(test)
    #print(test.index)
    a = np.dot(test,test.T) #dot product counts number of times cells are paired together
    #a.col_names =
    #print(a.shape)
    #print(a)
    for j,k in itertools.product(range(len(test.index)), range(len(test.index))):
        if a[j][k] == 0:
             continue
        cell1 = test.index[j]
        #print(cell1)
        cell2 = test.index[k]
        position_j_orig = original_list_indices.index(cell1)
        #print(position_j_orig)
        position_k_orig = original_list_indices.index(cell2)
        #print(position_k_orig)
        pairwise_coocurrence_val = pairwise_cooccurrence[position_j_orig][position_k_orig]
        pairwise_cooccurrence[position_j_orig][position_k_orig] = pairwise_coocurrence_val + a[j][k]
        #print(j)
        #print(k)

    print("completed run: " + str(i))

df = pd.DataFrame(pairwise_cooccurrence)
df.index = adata.obs_names
df.columns = adata.obs_names

df.to_csv("pairwise_clustering_results.csv")
