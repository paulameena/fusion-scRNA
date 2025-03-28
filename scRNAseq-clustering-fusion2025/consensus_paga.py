import scanpy as sc
import numpy as np
import pandas as pd
from sklearn.metrics import silhouette_score, silhouette_samples
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
import itertools
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import seaborn as sns

#reload adata
adata= sc.read_h5ad("adata_filtered_combined_feb2025.h5ad")
linkage_matrix = pd.read_csv("/mnt/pan/SOM_CCCC_JGS25/shultesp/scRNAseq-clustering-fusion2025/linkage.csv")

#Cut the dendrogram to get a certain number of clusters
cluster_labels = fcluster(linkage_matrix, t=12, criterion='maxclust')

#Step 1: add clusters to adata
adata.obs['consensus_cluster'] = cluster_labels
adata.obs['consensus_cluster']= adata.obs['consensus_cluster'].astype('category')


sc.tl.leiden( #leiden clustering algorithm is what determines number of clusters; tutorial has known number of cell types...
    adata,
    resolution=1,
    random_state=0,
    n_iterations=-1,
    directed=False,
)

#now create paga initialization and map clusters on umap and compare to original output...
#chatgpt helped me identify the correct packages to use

# Step 2: Compute the neighborhood graph
sc.pp.neighbors(adata)

# Step 3: Compute the PAGA graph
sc.tl.paga(adata, groups='consensus_cluster')  # Adjust to your cluster identification

# Step 4: Visualize the PAGA graph
sc.pl.paga(adata,save="_consensus_paga.png")

# Step 5: Run UMAP with PAGA initialization
sc.tl.umap(adata, init_pos='paga')  # Using the PAGA initialization

# Step 6: Visualize UMAP
sc.pl.umap(adata, color=['consensus_cluster', 'sample', 'leiden'], save="_consensus_umap.png")  # Plot UMAP with cluster colors

#compare to original leiden clusters

#also compare via diffusion maps
sc.pp.neighbors(adata)

# Step 3: Compute diffusion maps
sc.tl.diffmap(adata)

# Step 4: Visualize diffusion map
sc.pl.diffmap(adata, color=['consensus_cluster', 'leiden', 'sample'], save="_diffusion_figs.png")
